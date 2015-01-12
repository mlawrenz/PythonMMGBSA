import sys
import shutil
import amber_file_formatter
import subprocess
from subprocess import PIPE
import pickle
import optparse
import pylab
import numpy
import glob
import os
from os.path import split


# Helper Functions 

def get_pbbond_radii(model):
    if model==1:
        radii='mbondi'
    if model==2:
        radii='mbondi'
    if model==5:
        radii='mbondi2'
    if model==7:
        radii='mbondi'
    if model==8:
        radii='mbondi3'
    return radii

def run_linux_process(command):
    p=subprocess.Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()
    output, err=p.communicate()
    return output, err

    

def amber_mask_reducer(testmask, comparemask):
    # hack to workaround ambmask stack underflow if you specify every residue
    # also am keeping whole residues that any atom is within X of ligand
    newlist=[]
    residues_list=testmask.split(',')
    badlist=[n for (n, i) in enumerate(residues_list) if i=='']
    for i in badlist:
        residues_list.pop(i)
    residues_list=sorted([int(i) for i in residues_list])
    compareresidues_list=comparemask.split(',')
    badlist=[n for (n, i) in enumerate(compareresidues_list) if i=='']
    for i in badlist:
        compareresidues_list.pop(i)
    compareresidues_list=sorted([int(i) for i in compareresidues_list])
    # exclude whole moveable residues
    for res in compareresidues_list:
        if res in residues_list:
            index=residues_list.index(res)
            residues_list.pop(index)
    residue_list=numpy.array(sorted([int(i) for i in residues_list]))
    previous=residue_list[0]
    start=previous
    for res in residue_list[1:]:
        if res==(previous + 1):
            previous=res
            pass
        elif start-previous==0:
            newlist.append('%s' % previous)
            previous=res
            start=res
        else:
            newlist.append('%s-%s' % (start, previous))
            previous=res
            start=res
    if start-previous==0:
        newlist.append('%s' % start)
    else:
        newlist.append('%s-%s' % (start, previous))
    return ','.join(newlist)
        

def get_restraints(prot_radius, prmtop, inpcrd, ligrestraint=False):
    # select **whole** residues that are within prot_radius of the molecule
    # include MOL, with option of restraining it too
    # exclude explicit waters from all of this
    base_top=os.path.basename(prmtop) #workaround ambmask char limit
    base_crd=os.path.basename(inpcrd)
    wd=os.path.dirname(prmtop)
    origdir=os.getcwd()
    os.chdir(wd)
    if ligrestraint==True:
        # then set restraints to include MOL (allows just protein around
        # molecule to move)
        command="ambmask -p %s -c %s -find \"! :WAT & :MOL > @%s | :MOL\" | grep ATOM | awk '{print $5}' | grep -v \"\*\*\" | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
        converse="ambmask -p %s -c %s -find \"! :WAT & :MOL < @%s | ! :MOL\" | grep ATOM | awk '{print $5}' | grep -v \"\*\*\" | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
    else:
        # else do not include MOL, so it moves
        command="ambmask -p %s -c %s -find \"! :WAT & :MOL > @%s\" | grep ATOM |   grep -v \"\*\*\" | awk '{print $5}' | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
        converse="ambmask -p %s -c %s -find \"! :WAT & :MOL < @%s\" | grep ATOM |   grep -v \"\*\*\" | awk '{print $5}' | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
    mask=subprocess.check_output(command, shell=True)
    conversemask=subprocess.check_output(converse, shell=True)
    mask=amber_mask_reducer(mask, conversemask)
    name=base_top.split('.top')[0]
    # save restrained, and moveable residue atoms
    numpy.savetxt('%s_restraintresidues.txt' % name, [mask,], fmt='%s')
    numpy.savetxt('%s_moveableresidues.txt' % name, [conversemask,], fmt='%s')
    os.chdir(origdir)
    return mask
    

def get_simulation_commands(prefix, prmtop, inpcrd, outdir, gpu=False, restrain=False, nproc=16, mdrun=False):
    if gpu==True:
        print "USING GPU FOR MD, ASSUMING 4 DEVICES!"
        os.system('export CUDA_VISIBLE_DEVICES=0,1,2,3')
        program='pmemd.cuda'
    else:
        program='mpirun -n %s pmemd.MPI' % nproc
        print "running on %s processors" % nproc
    if restrain==True:
        if mdrun==True:
            command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -ref {4} -r {1}/{2}.rst -x {1}/{2}.mdcrd'.format(program, outdir, prefix, prmtop, inpcrd)
        else:
            command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -ref {4} -r {1}/{2}.rst'.format(program, outdir, prefix, prmtop, inpcrd)
    else:
        if mdrun==True:
            command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -r {1}/{2}.rst -x {1}/{2}.mdcrd'.format(program, outdir, prefix, prmtop, inpcrd)
        else:
            command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -r {1}/{2}.rst'.format(program, outdir, prefix, prmtop, inpcrd)
    return command


# Class for Building Amber Parametrized Molecule
class ambermol:
    '''sets up molecular parameters and input files for min (single point calc) or
MD, for processing with MMGB scores'''
    def __init__(self, jobname=None, protfile=None, ligfile=None, ligcharge=None, \
implicit=False, gbmodel=1, md=False, mdsteps=100000, mdseed=-1, maxcycles=50000, drms=0.1, nproc=8, gpu=False, \
prot_radius=None, restraint_k=10.0, ligrestraint=None):
        self.nproc=int(nproc)
        self.jobname=jobname
        self.gbmodel=int(gbmodel)
        self.radii=get_pbbond_radii(int(gbmodel))
        self.implicit=implicit
        self.maxcycles=int(maxcycles)
        self.drms=float(drms)
        self.md=md
        self.gpu=gpu
        self.mdsteps=mdsteps
        self.mdseed=mdseed
        self.protfile=os.path.abspath(protfile)
        self.ligfile=os.path.abspath(ligfile)
        self.ligand_name=os.path.basename(ligfile).split('.mol2')[0]
        print "--------------------------------------"
        print "SYSTEM SET UP-------------------------"
        # check ligand file for resname MOL
        command="more %s | awk '{if (NF==9) {print $8}}' | head -1" % self.ligfile
        output=subprocess.check_output(command, shell=True)
        output=output.rstrip('\n')
        if output!='MOL':
            print "NEED TO NAME LIGAND \"MOL\""
            sys.exit()
        # check for ligand charge
        if not ligcharge:
            print "CONFIRM LIGAND CHARGE"
            sys.exit()
        else:
            print "LIGAND CHARGE IS %s" % ligcharge
            self.ligcharge=int(ligcharge)
        self.charge_method='bcc'
        print "USING AM1-BCC CHARGE METHOD FOR LIGAND"
        print "USING MMGB=%s MODEL FOR MMGBSA CALC" % self.gbmodel
        # set tmp directories and main output directory self.gbdir
        self.antdir='%s/%s-antechamber-output' % (os.getcwd(), self.ligand_name)
        self.leapdir='%s/%s-leap-output' % (os.getcwd(), self.jobname)
        self.amberligfile='%s/%s.amber.mol2' % (self.antdir, self.ligand_name)
        if self.md==True:
            if implicit==True:
                self.mdprefix='gbmd-cpx'
                self.gbdir='%s/%s-implicit-gb%s-md' % (os.getcwd(), self.jobname, self.gbmodel)
            else:
                self.mdprefix='md-cpx'
                self.gbdir='%s/%s-explicit-gb%s-md' % (os.getcwd(), self.jobname, self.gbmodel)
        else:
            if implicit==True:
                self.gbdir='%s/%s-implicit-gb%s-min' % (os.getcwd(), self.jobname, self.gbmodel)
            else:
                self.gbdir='%s/%s-explicit-gb%s-min' % (os.getcwd(), self.jobname, self.gbmodel)
        if not os.path.exists(self.gbdir):
            os.mkdir(self.gbdir)
        #store minimized complex
        if self.implicit==True:
            self.mincpx='%s/implicit-cpx.rst' % self.gbdir
        else:
            self.mincpx='%s/min-cpx.rst' % self.gbdir
        if not os.environ['AMBERHOME']:
            print "AMBERHOME IS NOT SET"
            sys.exit()
        # set restraints
        if not prot_radius:
            print "NO PROTEIN RESTRAINTS USED"
            self.prot_radius=None
        else:
            self.prot_radius=prot_radius
            print "PROTEIN RESTRAINED AT RADIUS %s AROUND LIGAND" % self.prot_radius
        self.restraint_k=restraint_k
        print "RESTRAINT FORCE %s kcal/mol*A2" % restraint_k
        if not ligrestraint:
            self.ligrestraint=None
            print "NO LIGAND RESTRAINTS USED"
        else:
            self.ligrestraint=True
            print "RESTRAINING LIGAND ATOMS"
 
    def ptraj_rst_to_pdb(self, inconf, topo, dir):
        origdir=os.getcwd()
        os.chdir(dir)
        file=open('ptraj-pdb.in', 'w')
        file.write('trajin %s\n' % inconf)
        file.write('trajout %s.pdb pdb' % inconf.split('.rst')[0])
        file.close()
        command='cpptraj %s ptraj-pdb.in' % topo
        output, err=run_linux_process(command)
        prefix='ptraj'
        self.check_output(output, err, prefix, type='ptraj')
        return

    def check_output(self, output, err, prefix, type):
        types=['leap', 'ante', 'ptraj', 'md', 'MMGBSA']
        outdirs=[self.leapdir, self.antdir, self.gbdir, self.gbdir, self.gbdir]
        typedict=dict()
        for (t,odir) in zip(types, outdirs):
            typedict[t]=odir
        dir=typedict[type]
        errors=False
        if 'rror' in err or 'rror' in output:
            errors=True
            if type=='leap' and 'tl_getline()' in err.split('\n')[1]:
                if len(err.split('\n')) < 4:
                    errors=False
            else:
                logfile='%s/%s-%s-%s.err' % (dir, self.ligand_name, prefix, type)
                numpy.savetxt(logfile, output.split('\n'), fmt='%s')
        if 'Abort' in err or 'Abort' in output:
            errors=True
        logfile='%s/%s-%s-%s.err' % (dir, self.ligand_name, prefix, type)
        if len(err)!=0:
            numpy.savetxt(logfile, err.split('\n'), fmt='%s')
        if errors==True:
            print "ERRORS: for %s" % type
            print "check %s" % logfile
            sys.exit()
        elif 'Parameter file was not saved' in output:
            print "ERRORS: parameter file was not saved for %s" % type
            numpy.savetxt(logfile, output.split('\n'), fmt='%s')
            print "check %s" % logfile
            sys.exit()
        elif 'Could not open' in output:
            print "ERRORS:could not open file for %s" % type
            numpy.savetxt(logfile, output.split('\n'), fmt='%s')
            print "check %s" % logfile
            sys.exit()
        else:
            pass
        return

    def run_antechamber(self):
        #check to see if file exists
        if not os.path.exists(self.antdir):
            os.mkdir(self.antdir)
        origdir=os.getcwd()
        os.chdir(self.antdir)
        # run antechamber in antdir to avoid overwritten files by parallel jobs
        frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
        if os.path.exists(self.amberligfile) and os.path.exists(frcmodfile):
            print "CHARGED MOL2 AND FRCMOD FILES FOR %s in %s" % (self.ligand_name, self.antdir)
        else:
            print "--------------------------------------"
            print "RUNNING ANTECHAMBER ON LIGAND---------"
            print "OUTPUT %s-----------------------------" % self.antdir
            #make new mol2 file for ligand, including charges
            command='%s/bin/antechamber -i %s -fi mol2 -o %s -fo mol2 -c %s -nc %s' % (os.environ['AMBERHOME'], self.ligfile,  self.amberligfile, self.charge_method,self.ligcharge)
            output, err=run_linux_process(command)
            prefix='ante'
            self.check_output(output, err, prefix, type='ante')
            #make frcmod file to catch missing parameters
            prefix='parmchk'
            command='%s/bin/parmchk -i %s -o %s -f mol2' % (os.environ['AMBERHOME'], self.amberligfile, frcmodfile)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix, type='ante')
        
            os.system('rm ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF  rm sqm.*')
        os.chdir(origdir)
        return


    def run_leap(self):
        if not os.path.exists(self.leapdir):
            os.mkdir(self.leapdir)
        # run leap in antdir to avoid overwritten files by parallel jobs
        origdir=os.getcwd()
        os.chdir(self.leapdir)
        print "--------------------------------------"
        print "RUNNING LEAP FOR COMPLEX--------------"
        frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
        prefix='cpx'
        amber_file_formatter.write_leap(self.leapdir, prefix, self.ligand_name, self.radii, frcmodfile, self.amberligfile, self.protfile, complex=True, implicit=self.implicit)
        command='%s/bin/tleap -f %s/%s-%s-leaprc' % (os.environ['AMBERHOME'],
self.leapdir, self.ligand_name, prefix)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix, type='leap')
        os.chdir(origdir)
        return


    def simulation_guts(self, prefix, prmtop, inpcrd, mdrun=False):
        # write simulation run input files
        # pass in prefix, prmtop, and inpcrd appropriate for gb vs. explicit
        nproc=self.nproc
        if 'ligand' in prefix:
            restraint_atoms=None
            restrain=False
            if self.implicit==True:
                nproc=2
        else:
            restraint_atoms=self.restraint_atoms
            if self.restraint_atoms!=None:
                restrain=True
            else:
                restrain=False
        if mdrun==False:
            if self.implicit==True:
                print "--------------------------------------"
                print "RUNNING MINIMIZATION WITH IMPLICIT----"
                amber_file_formatter.write_simulation_input(md=False,dir=self.gbdir, prefix=prefix, implicit=self.implicit, gbmodel=self.gbmodel, restraint_k=self.restraint_k, restraint_atoms=restraint_atoms, maxcycles=self.maxcycles, drms=self.drms)
            else:
                print "RUNNING MINIMIZATION WITH EXPLICIT----"
                amber_file_formatter.write_simulation_input(md=False, dir=self.gbdir, prefix=prefix, restraint_atoms=restraint_atoms, restraint_k=self.restraint_k,maxcycles=self.maxcycles, drms=self.drms)
            command=get_simulation_commands(prefix, prmtop, inpcrd, self.gbdir, self.gpu, restrain, nproc)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix=prefix, type='md')
        else:
            print "--------------------------------------"
            if self.implicit==True:
                print "RUNNING MD SIMULATION WITH IMPLICIT----"
                amber_file_formatter.write_simulation_input(md=True, dir=self.gbdir, prefix=prefix,  implicit=self.implicit, gbmodel=self.gbmodel,\
restraint_atoms=restraint_atoms, restraint_k=self.restraint_k, steps=self.mdsteps, mdseed=self.mdseed)

            else:
                print "RUNNING MD SIMULATION WITH EXPLICIT---"
                amber_file_formatter.write_simulation_input(md=True, dir=self.gbdir, prefix=prefix,  restraint_atoms=restraint_atoms, \
restraint_k=self.restraint_k, steps=self.mdsteps, mdseed=self.mdseed)
            command=get_simulation_commands(prefix, prmtop, inpcrd, self.gbdir, self.gpu, restrain, nproc, mdrun=True)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='md', type='md')
        return
        

    def run_cpx_simulation(self):
        # set input variables and commands
        if self.implicit==True:
            prmtop='%s/%s-complex.top' % (self.leapdir, self.ligand_name)
            inpcrd='%s/%s-complex.crd' % (self.leapdir, self.ligand_name)
            prefix='implicit-cpx'
        else:
            prefix='min-cpx'
            prmtop='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
            inpcrd='%s/%s-complex.solv.crd' % (self.leapdir, self.ligand_name)
        if not os.path.exists(prmtop):
            print "MISSING TOPOLOGY FILE: %s" % prmtop 
            sys.exit()
        if not os.path.exists(inpcrd):
            print "MISSING COOR FILE: %s" % inpcrd
            sys.exit()
        if self.ligrestraint==True and self.prot_radius!=None:
            restraint_atoms=get_restraints(self.prot_radius, prmtop, inpcrd, ligrestraint=True)
        elif self.ligrestraint!=True and self.prot_radius!=None:
            restraint_atoms=get_restraints(self.prot_radius, prmtop, inpcrd, ligrestraint=False)
        else:
            restraint_atoms=None
        self.restraint_atoms=restraint_atoms
        self.simulation_guts(prefix, prmtop, inpcrd)
        print "self.mincpx is %s" % prefix
        if self.md==True:
            if not os.path.exists(self.mincpx):
                print "minimization failed, no %s" % self.mincpx
                sys.exit()
            inpcrd=self.mincpx
            self.simulation_guts(self.mdprefix, prmtop, inpcrd, mdrun=True)
        return

    def run_ligand_strain(self):
        # extract ligand from minimized complex structure with cpptraj
        print "--------------------------------------"
        print "SETTING UP LIGAND STRAIN CALC---------"
        if not os.path.exists(self.mincpx):
            print "minimization failed, no %s" % self.mincpx
            sys.exit()
        base=os.path.basename(self.mincpx)
        filename='%s/getligand.ptraj' % self.gbdir
        if self.implicit==True:
            minligand='%s/ligand_in_cpx.rst' % self.gbdir #rst for sim start
            # ptraj file iteself in same dir as inconf and outconf
            # but calling it from one above
            amber_file_formatter.write_ptraj_strip(filename, self.mincpx,minligand)
            #extract ligand to get restart file directly for MMPBSA
            command='cpptraj {0}/{1}-complex.top {2}'.format(self.leapdir, self.ligand_name, filename)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='strip-ligand', type='ptraj')
            #for running simulation
            prefix='implicit-ligand'
            inpcrd=minligand
            prmtop='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
            self.ptraj_rst_to_pdb(minligand, prmtop, self.gbdir)
        else:
            # rebuild ligand topology only if using explicit solvent
            # has to load in a mol2
            minligand='%s/ligand_in_cpx.mol2' % self.gbdir
            amber_file_formatter.write_ptraj_strip(filename, self.mincpx, minligand)
            command='cpptraj {0}/{1}-complex.solv.top {2}'.format(self.leapdir, self.ligand_name, filename)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='strip-ligand', type='ptraj')
            # rebuild ligand topology 
            frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
            leap_prefix='ligresolv'
            amber_file_formatter.write_leap(dir=self.leapdir, prefix=leap_prefix, ligand_name=self.ligand_name, radii=self.radii, frcmodfile=frcmodfile, newligfile=minligand, complex=False, implicit=False)
            command='{0}/bin/tleap -f {1}/{2}-{3}-leaprc'.format(os.environ['AMBERHOME'], self.leapdir, self.ligand_name, leap_prefix)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix=leap_prefix, type='leap')
            print "RAN LEAP for RESOLVATING LIGAND"
            os.system('mv %s/%s-ligand.solv.crd %s/ligand_in_cpx.crd' % (self.leapdir, self.ligand_name, self.gbdir))
            #for running simulation
            prefix='min-ligand'
            inpcrd='%s/ligand_in_cpx.crd' % self.gbdir
            prmtop='%s/%s-ligand.solv.top' % (self.leapdir, self.ligand_name)
        # minimize ligand in solution (or with GB solvent)
        # see diff in complex energy and solvated energy
        self.simulation_guts(prefix, prmtop, inpcrd)
        print "MINIMIZED LIGAND"
        self.ptraj_rst_to_pdb('%s.rst' % prefix, prmtop, self.gbdir)
        self.run_mmgbsa(complex=False)
        print "MMGB CALC FINISHED ON LIGAND"
        return

    def mmgbsa_guts(self, prefix, start, finish, solvcomplex, complex, traj, interval=1, protein=None, ligand=None):
        # should be in gbdir here
        inputfile='%s-mmgb.in' % prefix
        # use reduced number of processes so don't run out of memory
        pb_processes=round(self.nproc/2.0)
        if self.md==True:   
            program='mpirun -n {0} {1}/bin/MMPBSA.py.MPI'.format(pb_processes, os.environ['AMBERHOME'])
        else:
            program='{0}/bin/MMPBSA.py'.format(os.environ['AMBERHOME'])
        amber_file_formatter.write_mmgbsa_input(inputfile, self.gbmodel, start, interval, finish)
        # use MMGBSA.py in Amber14 to run MMGB free energy difference calcs for complex
        if protein!=None and ligand!=None:
            print "--------------------------------------"
            print "RUNNING COMPLEX MMGBSA CALC-----------"
            if self.implicit==True:
                command='{0} -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -cp {4} \
                 -rp {5} -lp {6} -y {7}'.format(program, inputfile, self.ligand_name, prefix, complex, protein, ligand, traj)
            else:
                command='{0} -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -sp {4} \
                -cp {5} -rp {6} -lp {7} -y {8}'.format(program, inputfile, self.ligand_name, prefix, solvcomplex, complex, protein, ligand, traj)
        else:
            print "--------------------------------------"
            print "RUNNING LIGAND MMGBSA CALC------------"
            program='{0}/bin/MMPBSA.py'.format(os.environ['AMBERHOME'])
            if self.implicit==True:
                command='{0} -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -cp {4} -y \
{5}'.format(program, inputfile, self.ligand_name, prefix, complex, traj)
            else:
                command='{0} -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -sp {4} \
                -cp {5} -y {6}'.format(program, inputfile, self.ligand_name, prefix, solvcomplex, complex, traj)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix, type='MMGBSA')
        return

    def run_mmgbsa(self, complex=True):
        # setup files for MMGBSA.py in Amber14 to run MMGB free energy difference calcs for complex
        origdir=os.getcwd()
        # make sure do not mixup extraneous files with parallel jobs
        os.chdir(self.gbdir)
        if complex==True:
            prefix='cpx'
            start=0
            interval=1
            if self.md==True:
                finish=1000000000 # MMPBSA.py will reduce to total frames
                if self.implicit==True:
                    traj='gbmd-cpx.mdcrd' 
                    solvcomplex=None
                else:
                    traj='md-cpx.mdcrd' 
                    solvcomplex='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
            else:
                finish=1
                if self.implicit==True:
                    traj='implicit-cpx.rst' 
                    solvcomplex=None
                else:
                    traj='min-cpx.rst' 
                    solvcomplex='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
            complex='%s/%s-complex.top' % (self.leapdir, self.ligand_name)
            protein='%s/%s-protein.top' % (self.leapdir, self.ligand_name)
            ligand='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, traj, interval=1, protein=protein, ligand=ligand)
        else:
            #setup files for MMGBSA.py in Amber14 to run MMGB free energy for ligand strain
            start=0
            interval=1
            finish=1
            if self.implicit==True:
                solvcomplex=None
                complex='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
                initial_traj='ligand_in_cpx.rst' 
                final_traj='implicit-ligand.rst' 
            else:
                solvcomplex='%s/%s-ligand.solv.top' % (self.leapdir, self.ligand_name)
                complex='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
                initial_traj='ligand_in_cpx.crd' 
                final_traj='min-ligand.rst' 
            # first get initial GB energy of ligand in complex
            prefix='ligcpx'
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, initial_traj, interval=1)
            # next GB energy of ligand minimized in solution
            prefix='ligsolv'
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, final_traj, interval=1)
        os.chdir(origdir)
        return

    def clean(self):
        os.system('cp %s/*pdb %s' % (self.leapdir, self.gbdir))
        os.system('cp %s/*-complex.solv* %s' % (self.leapdir, self.gbdir))
        os.system('cp %s/*-ligand.solv* %s' % (self.leapdir, self.gbdir))
        os.system('cp %s/*.amber.mol2 %s/' % (self.antdir, self.gbdir))
        # convert all restart files into PDB files
        if self.implicit==True:
            prmtop='%s/%s-complex.top' % (self.leapdir, self.ligand_name)
        else:
            prmtop='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
        if os.path.exists(self.mincpx):
            self.ptraj_rst_to_pdb(self.mincpx, prmtop, self.gbdir)
        else:
            print "MISSING MINIMIZED COMPLEX STRUCTURE: ", self.mincpx
            sys.exit()
        if self.md==True:
            if os.path.exists('%s.rst' % self.mdprefix):
                self.ptraj_rst_to_pdb('%s.rst' % self.mdprefix, prmtop, self.gbdir)
            else:
                print "MISSING MD COMPLEX STRUCTURE: ", self.mincpx
                sys.exit()
        shutil.rmtree(self.antdir)
        shutil.rmtree(self.leapdir)
        os.system('rm %s/*in' % self.gbdir)
        os.system('rm %s/*rst' % self.gbdir)
        os.system('rm %s/*out' % self.gbdir)
        os.system('rm %s/*ptraj*' % self.gbdir)
        print "FINISHED CLEANING: use -debug if you want to check intermediate files"
        

    def print_table(self):
        dir=self.gbdir
        all_errors=['MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB', 'E_surf', 'E_lig']
        all_values=dict()
        all_errors=dict()
        files=glob.glob('%s/*-cpx-*FINAL*' % dir)
        if len(files) ==0:
            print "MISSING MMPBSA OUTPUT"
            sys.exit()
        allcomponents=['']
        for file in files:
            base=os.path.basename(file)
            ligand=base.split('-cpx-')[0]
            if ligand not in all_values.keys():
                all_values[ligand]=dict()
                all_errors[ligand]=dict()
            collect=False
            print "on %s" % file
            fhandle=open(file)
            for line in fhandle.readlines():
                if 'Differences' in line:
                    collect=True
                if collect==True:
                    if 'DELTA TOTAL' in line:
                        all_values[ligand]['MMGB']=float(line.split()[2])
                        all_errors[ligand]['MMGB']=float(line.split()[3])
                    if 'VDWAALS' in line:
                        all_values[ligand]['vdW']=float(line.split()[1])
                        all_errors[ligand]['vdW']=float(line.split()[2])
                    if 'EEL' in line and '1-4 EEL' not in line:
                        all_values[ligand]['eel_inter']=float(line.split()[1])
                        all_errors[ligand]['eel_inter']=float(line.split()[2])
                    if 'EGB' in line:
                        all_values[ligand]['EGB']=float(line.split()[1])
                        all_errors[ligand]['EGB']=float(line.split()[2])
                    if 'ESURF' in line:
                        all_values[ligand]['E_surf']=float(line.split()[1])
                        all_errors[ligand]['E_surf']=float(line.split()[2])
            all_values[ligand]['eel/EGB']=all_values[ligand]['eel_inter']+all_values[ligand]['EGB']
            all_errors[ligand]['eel/EGB']=numpy.sqrt(all_errors[ligand]['eel_inter']**2 + all_errors[ligand]['EGB']**2)
        states=['ligcpx', 'ligsolv'] # cpx - solv = strain
        components=dict()
        for ligandstate in states:
            files=glob.glob('%s/*-%s-*FINAL*' % (dir, ligandstate))
            if len(files) ==0:
                print "MISSING LIGAND STRAIN CALC"
                sys.exit()
            values=[]
            errors=[]
            components[ligandstate]=dict()
            for file in files:
                base=os.path.basename(file)
                ligand=base.split('-%s-' % ligandstate)[0]
                components[ligandstate][ligand]=dict()
                fhandle=open(file)
                for line in fhandle.readlines():
                    if 'TOTAL' in line:
                        components[ligandstate][ligand]['value']=float(line.split()[1])
                        components[ligandstate][ligand]['err']=float(line.split()[2])
                    if ligandstate=='ligsolv' and 'G gas' in line:
                        all_values[ligand]['E_lig']=line.split()[2]
                        all_errors[ligand]['E_lig']=line.split()[3]
        for ligand in components['ligcpx'].keys():
            all_values[ligand]['strain']=components['ligcpx'][ligand]['value']-components['ligsolv'][ligand]['value']
            all_errors[ligand]['strain']=error=numpy.sqrt(components['ligcpx'][ligand]['err']**2+components['ligsolv'][ligand]['err']**2)
        ligands=all_values.keys()
        sorted_ligands=sorted(ligands, key=lambda x: all_values[x]['MMGB'])
        ohandle=open('%s/results.tbl' % dir, 'w')
        #format with 9 max column width
        entry='{0:<9} {1:<9} {2:<9} {3:<9} {4:<9} {5:<9} {6:<9} {7:<9} {8:<9} {9:<9}'.format('root', 'MMGB+str', 'MMGB', \
'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB' , 'E_surf',  'E_lig')
        entry=''.join([ entry, '\n'])
        ohandle.write(entry)
        print entry
        keyorder=['MMGB+str', 'MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB' , 'E_surf',  'E_lig'] 
        for ligand in sorted_ligands:
            all_values[ligand]['MMGB+str']=all_values[ligand]['MMGB']+all_values[ligand]['strain']
            name='%-10s' % ligand
            entry=''.join(['%-10.2f' % round(float(all_values[ligand][x]), 2) for x in keyorder])
            entry=''.join([name, entry, '\n'])
            ohandle.write(entry)
            print entry


