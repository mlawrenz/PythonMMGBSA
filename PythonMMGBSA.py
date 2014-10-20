import sys
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

def amber_mask_reducer(testmask):
    # hack to workaround ambmask stack underflow if you specify every residue
    residues_list=testmask.split(',')
    newlist=[]
    badlist=[n for (n, i) in enumerate(residues_list) if i=='']
    for i in badlist:
        residues_list.pop(i)
    residue_list=numpy.array(sorted([int(i) for i in residues_list]))
    previous=residue_list[0]
    start=previous
    for res in residue_list[1:]:
        if res==(previous + 1):
            previous=res
            pass
        else:
            newlist.append('%s-%s' % (start, previous))
            previous=res
            start=res
    newlist.append('%s-%s' % (start, previous))
    return ','.join(newlist)
        

def get_restraints(prot_radius, prmtop, inpcrd, ligrestraint=False):
    # select all atoms in residues that are within prot_radius of the molecule
    # include MOL, with option of restraining it too
    base_top=os.path.basename(prmtop) #workaround ambmask char limit
    base_crd=os.path.basename(inpcrd)
    wd=os.path.dirname(prmtop)
    origdir=os.getcwd()
    os.chdir(wd)
    if ligrestraint==True:
        # then set restraints to include MOL (allows just protein around
        # molecule to move)
        command="ambmask -p %s -c %s -find \":MOL > @%s | :MOL\" | grep ATOM | awk '{print $5}' | grep -v \"\*\*\" | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
        converse="ambmask -p %s -c %s -find \":MOL < @%s | ! :MOL\" | grep ATOM | awk '{print $5}' | grep -v \"\*\*\" | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
    else:
        # else do not include MOL, so it moves
        command="ambmask -p %s -c %s -find \":MOL > @%s\" | grep ATOM |   grep -v \"\*\*\" | awk '{print $5}' | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
        converse="ambmask -p %s -c %s -find \":MOL < @%s\" | grep ATOM |   grep -v \"\*\*\" | awk '{print $5}' | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, prot_radius)
    mask=subprocess.check_output(command, shell=True)
    conversemask=subprocess.check_output(converse, shell=True)
    mask=amber_mask_reducer(mask)
    conversemask=amber_mask_reducer(conversemask)
    #save restrained residue atoms
    command="ambmask -p %s -c %s -find :%s" % (base_top, base_crd, mask)
    output=subprocess.check_output(command, shell=True)
    name=prmtop.split('.top')[0]
    numpy.savetxt('%s_restraintresidues.txt' % name, output.split('\n'), fmt='%s')
    #save moveable residue atoms
    command="ambmask -p %s -c %s -find :%s" % (base_top, base_crd, conversemask)
    output=subprocess.check_output(command, shell=True)
    numpy.savetxt('%s_moveableresidues.txt' % name, output.split('\n'), fmt='%s')
    os.chdir(origdir)
    return mask
    
        
# Class for Building Amber Parametrized Molecule
class ambermol:
    '''sets up molecular parameters and input files for min (single point calc) or
MD, for processing with MMGB scores'''
    def __init__(self, jobname, protfile=None, ligfile=None, prot_radius=None, ligrestraint=None, charge_method=None, ligcharge=None, gbmin=False, gbmodel=5, restraint_k=5.0, md=False, mdsteps=50000, maxcycles=50000, drms=0.1, gpu=False, verbose=False):
        self.jobname=jobname
        self.verbose=verbose
        self.restraint_k=restraint_k
        print "RESTRAINT FORCE %s kcal/mol*A2" % restraint_k
        self.gbmodel=int(gbmodel)
        self.radii=get_pbbond_radii(int(gbmodel))
        self.gbmin=gbmin
        print "--------------------------------------"
        print "SYSTEM SET UP-------------------------"
        print "USING MMGB=%s MODEL" % self.gbmodel
        self.protfile='%s/%s' % (os.getcwd(), protfile)
        self.ligfile='%s/%s' % (os.getcwd(), ligfile)
        command="more %s | awk '{if (NF==9) {print $8}}' | head -1" % self.ligfile
        output=subprocess.check_output(command, shell=True)
        output=output.rstrip('\n')
        if output!='MOL':
            print "NEED TO NAME LIGAND \"MOL\""
            sys.exit()
        self.ligand_name=os.path.basename(ligfile).split('.mol2')[0]
        self.antdir='%s/%s-antechamber-output' % (os.getcwd(), self.jobname)
        self.leapdir='%s/%s-leap-output' % (os.getcwd(), self.jobname)
        self.amberligfile='%s/%s.amber.mol2' % (self.antdir, self.ligand_name)
        self.md=md
        self.gpu=gpu
        self.mdsteps=mdsteps
        if self.md==True:
            self.gbdir='%s-mmgb%s-md' % (self.jobname, self.gbmodel)
        else:
            self.gbdir='%s-mmgb%s-min' % (self.jobname, self.gbmodel)
        if not os.path.exists(self.gbdir):
            os.mkdir(self.gbdir)
        self.maxcycles=int(maxcycles)
        self.drms=float(drms)
        if not os.environ['AMBERHOME']:
            print "AMBERHOME IS NOT SET"
            sys.exit()
        if not prot_radius:
            print "NO PROTEIN RESTRAINTS USED"
            self.prot_radius=None
        else:
            self.prot_radius=prot_radius
            print "PROTEIN RESTRAINED AT RADIUS %s AROUND LIGAND" % self.prot_radius
        if not ligrestraint:
            self.ligrestraint=None
            print "NO LIGAND RESTRAINTS USED"
        else:
            self.ligrestraint=True
            print "RESTRAINING LIGAND ATOMS"
        if not charge_method:
            print "USING AM1-BCC CHARGE METHOD"
            self.charge_method='bcc'
        else: 
            self.charge_method=charge_method
        if not ligcharge:
            print "CONFIRM LIGAND CHARGE"
            sys.exit()
        else:
            print "LIGAND CHARGE IS %s" % ligcharge
            self.ligcharge=int(ligcharge)
 
    def get_simulation_commands(self, prefix, prmtop, inpcrd, restrain=False, nproc=16, mdrun=False):
        if self.gpu==True:
            print "USING GPU FOR MD"
            os.system('export CUDA_VISIBLE_DEVICES=0')
            program='pmemd.cuda'
        else:
            program='mpirun -n %s pmemd.MPI' % nproc
        if restrain==True:
            if mdrun==True:
                command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -ref {4} -r {1}/{2}.rst -x {1}/{2}.mdcrd'.format(program, self.gbdir, prefix, prmtop, inpcrd)
            else:
                command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -ref {4} -r {1}/{2}.rst'.format(program, self.gbdir, prefix, prmtop, inpcrd)
        else:
            if mdrun==True:
                command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -r {1}/{2}.rst -x {1}/{2}.mdcrd'.format(program, self.gbdir, prefix, prmtop, inpcrd)
            else:
                command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -r {1}/{2}.rst'.format(program, self.gbdir, prefix, prmtop, inpcrd)
        return command

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
        if self.verbose==True:
            logfile='%s/%s-%s-%s.log' % (dir, self.ligand_name, prefix, type)
            numpy.savetxt(logfile, output.split('\n'), fmt='%s')
            print "OUTPUT SAVED IN %s" % logfile
        logfile='%s/%s-%s-%s.err' % (dir, self.ligand_name, prefix, type)
        if len(err)!=0:
            numpy.savetxt(logfile, err.split('\n'), fmt='%s')
        if errors==True:
            print "ERRORS: for %s" % type
            print "check %s" % logfile
            sys.exit()
        elif 'Parameter file was not saved' in output:
            print "ERRORS: for %s" % type
            print "check %s" % logfile
            sys.exit()
        elif 'Could not open' in output:
            print "ERRORS: for %s" % type
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
            print "ALREADY HAVE CHARGED MOL2 AND FRCMOD FILES FOR %s" % self.ligand_name
        else:
            print "--------------------------------------"
            print "RUNNING ANTECHAMBER ON LIGAND---------"
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
        amber_file_formatter.write_leap(self.leapdir, prefix, self.ligand_name, self.radii, frcmodfile, self.amberligfile, self.protfile, complex=True, gbmin=self.gbmin)
        command='%s/bin/tleap -f %s/%s-%s-leaprc' % (os.environ['AMBERHOME'],
self.leapdir, self.ligand_name, prefix)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix, type='leap')
        os.chdir(origdir)
        return


    def simulation_guts(self, prefix, prmtop, inpcrd, mdrun=False):
        # write simulation run input files
        # pass in prefix, prmtop, and inpcrd appropriate for gb vs. explicit
        nproc=8
        if 'ligand' in prefix:
            restraint_atoms=None
            restrain=False
            if self.gbmin==True:
                nproc=2
        else:
            restraint_atoms=self.restraint_atoms
            if self.restraint_atoms!=None:
                restrain=True
            else:
                restrain=False
        if mdrun==False:
            if self.gbmin==True:
                print "--------------------------------------"
                print "RUNNING MINIMIZATION WITH IMPLICIT----"
                amber_file_formatter.write_simulation_input(md=False,dir=self.gbdir, prefix=prefix, gbmin=self.gbmin, gbmodel=self.gbmodel, restraint_k=self.restraint_k, restraint_atoms=restraint_atoms, maxcycles=self.maxcycles, drms=self.drms)
            else:
                print "RUNNING MINIMIZATION WITH EXPLICIT----"
                amber_file_formatter.write_simulation_input(md=False, dir=self.gbdir, prefix=prefix, restraint_atoms=restraint_atoms, restraint_k=self.restraint_k,maxcycles=self.maxcycles, drms=self.drms)
            command=self.get_simulation_commands(prefix, prmtop, inpcrd, restrain, nproc)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix=prefix, type='md')
        else:
            print "--------------------------------------"
            if self.gbmin==True:
                print "RUNNING MD SIMULATION WITH IMPLICIT----"
                amber_file_formatter.write_simulation_input(md=True, dir=self.gbdir, prefix=prefix,  gbmin=self.gbmin, gbmodel=self.gbmodel, restraint_atoms=restraint_atoms, restraint_k=self.restraint_k, steps=self.mdsteps)
            else:
                print "RUNNING MD SIMULATION WITH EXPLICIT---"
                amber_file_formatter.write_simulation_input(md=True, dir=self.gbdir, prefix=prefix,  restraint_atoms=restraint_atoms, restraint_k=self.restraint_k, steps=self.mdsteps)
            command=self.get_simulation_commands(prefix, prmtop, inpcrd, restrain, nproc, mdrun=True)

            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='md', type='md')
        return
        

    def run_cpx_simulation(self):
        # set input variables and commands
        if self.gbmin==True:
            prmtop='%s/%s-complex.top' % (self.leapdir, self.ligand_name)
            inpcrd='%s/%s-complex.crd' % (self.leapdir, self.ligand_name)
            prefix='gbmin-cpx'
            if self.md==True:
                mdprefix='gbmd-cpx'
        else:
            prefix='min-cpx'
            prmtop='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
            inpcrd='%s/%s-complex.solv.crd' % (self.leapdir, self.ligand_name)
            if self.md==True:
                mdprefix='md-cpx'
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
        #store minimized complex
        self.mincpx='%s/%s.rst' % (self.gbdir, prefix)
        print "self.mincpx is %s" % prefix
        if self.md==True:
            inpcrd=self.mincpx
            self.simulation_guts(mdprefix, prmtop, inpcrd, mdrun=True)
        return

    def run_ligand_strain(self):
        # extract ligand from minimized complex structure with cpptraj
        print "--------------------------------------"
        print "SETTING UP LIGAND STRAIN CALC---------"
        base=os.path.basename(self.mincpx)
        filename='%s/getligand.ptraj' % self.gbdir
        if self.gbmin==True:
            minligand='%s/ligandonly.rst' % self.gbdir #rst for sim start
            # ptraj file iteself in same dir as inconf and outconf
            # but calling it from one above
            amber_file_formatter.write_ptraj_strip(filename, self.mincpx,minligand)
            #extract ligand to get restart file directly for MMPBSA
            command='cpptraj {0}/{1}-complex.top {2}'.format(self.leapdir, self.ligand_name, filename)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='strip-ligand', type='ptraj')
            #for running simulation
            prefix='gbmin-ligand'
            inpcrd=minligand
            prmtop='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
        else:
            # rebuild ligand topology only if using explicit solvent
            # has to load in a mol2
            minligand='%s/ligandonly.mol2' % self.gbdir
            amber_file_formatter.write_ptraj_strip(filename, self.mincpx, minligand)
            command='cpptraj {0}/{1}-complex.solv.top {2}'.format(self.leapdir, self.ligand_name, filename)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='strip-ligand', type='ptraj')
            # rebuild ligand topology 
            frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
            leap_prefix='ligresolv'
            amber_file_formatter.write_leap(dir=self.leapdir, prefix=leap_prefix, ligand_name=self.ligand_name, radii=self.radii, frcmodfile=frcmodfile, newligfile=minligand, complex=False, gbmin=False)
            command='{0}/bin/tleap -f {1}/{2}-{3}-leaprc'.format(os.environ['AMBERHOME'], self.leapdir, self.ligand_name, leap_prefix)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix=leap_prefix, type='leap')
            print "RAN LEAP for RESOLVATING LIGAND"
            os.system('mv %s/%s-ligand.solv.crd %s/ligandonly.crd' % (self.leapdir, self.ligand_name, self.gbdir))
            #for running simulation
            prefix='min-ligand'
            inpcrd='%s/ligandonly.crd' % self.gbdir
            prmtop='%s/%s-ligand.solv.top' % (self.leapdir, self.ligand_name)
        # minimize ligand in solution (or with GB solvent)
        # see diff in complex energy and solvated energy
        self.simulation_guts(prefix, prmtop, inpcrd)
        print "MINIMIZED LIGAND"
        self.run_mmgbsa(complex=False)
        print "MMGB CALC FINISHED ON LIGAND"
        return

    def mmgbsa_guts(self, prefix, start, finish, solvcomplex, complex, traj, interval=1, protein=None, ligand=None):
        # should be in gbdir here
        inputfile='%s-mmgb.in' % prefix
        amber_file_formatter.write_mmgbsa_input(inputfile, self.gbmodel, start, interval, finish)
        # use MMGBSA.py in Amber14 to run MMGB free energy difference calcs for complex
        if protein!=None and ligand!=None:
            print "--------------------------------------"
            print "RUNNING COMPLEX MMGBSA CALC-----------"
            if self.gbmin==True:
                command='{0}/bin/MMPBSA.py -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -cp {4} -rp {5} -lp {6} -y {7}'.format(os.environ['AMBERHOME'], inputfile, self.ligand_name, prefix, complex, protein, ligand, traj)
            else:
                command='{0}/bin/MMPBSA.py -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -sp {4} -cp {5} -rp {6} -lp {7} -y {8}'.format(os.environ['AMBERHOME'], inputfile, self.ligand_name, prefix, solvcomplex, complex, protein, ligand, traj)
        else:
            print "--------------------------------------"
            print "RUNNING LIGAND MMGBSA CALC------------"
            if self.gbmin==True:
                command='{0}/bin/MMPBSA.py -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -cp {4} -y {5}'.format(os.environ['AMBERHOME'], inputfile, self.ligand_name, prefix, complex, traj)
            else:
                command='{0}/bin/MMPBSA.py -i {1} -o {2}-{3}-FINAL_MMPBSA.dat -sp {4} -cp {5} -y {6}'.format(os.environ['AMBERHOME'], inputfile, self.ligand_name, prefix, solvcomplex, complex, traj)
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
                if self.gbmin==True:
                    traj='gbmd-cpx.mdcrd' 
                    solvcomplex=None
                else:
                    traj='md-cpx.mdcrd' 
                    solvcomplex='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
            else:
                finish=1
                if self.gbmin==True:
                    traj='gbmin-cpx.rst' 
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
            if self.gbmin==True:
                solvcomplex=None
                complex='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
                initial_traj='ligandonly.rst' 
                final_traj='gbmin-ligand.rst' 
            else:
                solvcomplex='%s/%s-ligand.solv.top' % (self.leapdir, self.ligand_name)
                complex='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
                initial_traj='ligandonly.crd' 
                final_traj='min-ligand.rst' 
            # first get initial GB energy of ligand in complex
            prefix='ligcpx'
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, initial_traj, interval=1)
            # next GB energy of ligand minimized in solution
            prefix='ligsolv'
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, final_traj, interval=1)
        os.chdir(origdir)
        return

    def print_table(self):
        dir=self.gbdir
        all_errors=['MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB', 'E_surf', 'E_lig']
        all_values=dict()
        all_errors=dict()
        files=glob.glob('%s/*-cpx-*FINAL*' % dir)
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
        ohandle=open('%s/sorted_results.tbl' % dir, 'w')
        formatkeyorder=['name  MMGB+str', 'MMGB  strain  vdW  eel_inter  eel/EGB  EGB  E_surf  E_lig'] 
        entry=''.join(['%s\t' % x for x in formatkeyorder])
        entry=''.join([ entry, '\n'])
        ohandle.write(entry)
        print entry
        keyorder=['MMGB+str', 'MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB' , 'E_surf',  'E_lig'] 
        for ligand in sorted_ligands:
            all_values[ligand]['MMGB+str']=all_values[ligand]['MMGB']+all_values[ligand]['strain']
            name='%s\t\t' % ligand
            entry=''.join(['%0.2f\t' % round(float(all_values[ligand][x]), 2) for x in keyorder])
            entry=''.join([name, entry, '\n'])
            ohandle.write(entry)
            print entry


