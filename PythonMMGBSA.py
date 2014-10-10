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

def get_belly_restraints(protein_radius, prmtop, inpcrd, ligand_restraints=False):
    #select all atoms in residues that are within protein_radius of the molecule
    # include MOL, with option of restraining it too
    base_top=os.path.basename(prmtop) #workaround ambmask char limit
    base_crd=os.path.basename(inpcrd)
    wd=os.path.dirname(prmtop)
    origdir=os.getcwd()
    os.chdir(wd)
    if ligand_restraints==True:
        # then set belly to not include MOL (relaxes protein around molecule)
        command="ambmask -p %s -c %s -find \":MOL < @%s & ! :MOL\" | grep ATOM | awk '{print $5}' | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, protein_radius)
    else:
        # else set belly to include MOL, so it moves
        command="ambmask -p %s -c %s -find \":MOL < @%s\" | grep ATOM | awk '{print $5}' | sort | uniq | tr \"\n\" \", \"" % (base_top, base_crd, protein_radius)
    mask=subprocess.check_output(command, shell=True)
    command="ambmask -p %s -c %s -find :%s" % (base_top, base_crd, mask)
    output=subprocess.check_output(command, shell=True)
    #save belly residue atoms
    name=prmtop.split('.top')[0]
    numpy.savetxt('%s_bellyresidues.txt' % name, output.split('\n'), fmt='%s')
    os.chdir(origdir)
    return mask
    
        
# Class for Building Amber Parametrized Molecule
class ambermol:
    '''sets up molecular parameters and input files for min (single point calc) or
MD, for processing with MMGB scores'''
    def __init__(self, proteinfile=None, ligandfile=None, protein_radius=None, ligand_restraints=None, charge_method=None, ligand_charge=None, gbmin=False, gb_model=5, restraint_k=0.5, md=False, md_steps=100000, maxcycles=50000, gpu=False, verbose=False):
        self.verbose=verbose
        self.restraint_k=restraint_k
        self.gb_model=int(gb_model)
        self.radii=get_pbbond_radii(int(gb_model))
        self.gbmin=gbmin
        print "--------------------------------------"
        print "SYSTEM SET UP-------------------------"
        print "USING MMGB=%s MODEL" % self.gb_model
        self.proteinfile='%s/%s' % (os.getcwd(), proteinfile)
        self.ligandfile='%s/%s' % (os.getcwd(), ligandfile)
        command="more %s | awk '{if (NF==9) {print $8}}' | head -1" % self.ligandfile
        output=subprocess.check_output(command, shell=True)
        output=output.rstrip('\n')
        if output!='MOL':
            print "NEED TO NAME LIGAND \"MOL\""
            sys.exit()
        self.ligand_name=os.path.basename(ligandfile).split('.mol2')[0]
        self.antdir='%s/antechamber-output' % os.getcwd()
        self.leapdir='%s/leap-output' % os.getcwd()
        self.amberligandfile='%s/%s.amber.mol2' % (self.antdir, self.ligand_name)
        self.md=md
        if self.md==True:
            self.gbdir='mmgb%s-md' % self.gb_model
        else:
            self.gbdir='mmgb%s-min' % self.gb_model
        if not os.path.exists(self.gbdir):
            os.mkdir(self.gbdir)
        if gpu==True:
            os.system('export CUDA_VISIBLE_DEVICES=0')
            self.program='pmemd.cuda'
        else:
            self.program='mpirun -n 8 pmemd.MPI'
        self.maxcycles=int(maxcycles)
        if not os.environ['AMBERHOME']:
            print "AMBERHOME IS NOT SET"
            sys.exit()
        if not protein_radius:
            print "NO PROTEIN RESTRAINTS USED"
            self.protein_radius=None
        else:
            self.protein_radius=protein_radius
            print "PROTEIN RESTRAINED AT RADIUS %s AROUND LIGAND" % self.protein_radius
            print "RESTRAINT FORCE %s kcal/mol*A2" % restraint_k
        if not ligand_restraints:
            self.ligand_restraints=None
            print "NO LIGAND RESTRAINTS USED"
        else:
            self.ligand_restraints=True
            print "RESTRAINING LIGAND ATOMS"
            print "RESTRAINT FORCE %s kcal/mol*A2" % restraint_k
        if not charge_method:
            print "USING AM1-BCC CHARGE METHOD"
            self.charge_method='bcc'
        else: 
            self.charge_method=charge_method
        if not ligand_charge:
            # to calc change with umt:
            # cat conf_0001.mol2 | umt .mol2 .mol2 > tmp.mol2
            # command="grep \"\<0\>\" tmp.mol2  | sed 1d | awk '{sum+=} END {print sum}'"
            #output=subprocess.self.check_output(command, shell=True)
            # rm tmp.mol2
            print "ASSUMING LIGAND CHARGE IS ZERO"
            self.ligand_charge=0
        else:
            print "LIGAND CHARGE IS %s" % ligand_charge
            self.ligand_charge=int(ligand_charge)
 
    def get_simulation_commands(self, prefix, prmtop, inpcrd, restrain=False):
        if restrain==True:
            command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -ref {4} -r {1}/{2}.rst'.format(self.program, self.gbdir, prefix, prmtop, inpcrd)
        else:
            command='{0} -O -i {1}/{2}.in -o {1}/{2}.out -p {3} -c {4} -r {1}/{2}.rst'.format(self.program, self.gbdir, prefix, prmtop, inpcrd)
        return command

    def run_antechamber(self):
        #check to see if file exists
        if not os.path.exists(self.antdir):
            os.mkdir(self.antdir)
        print "--------------------------------------"
        print "RUNNING ANTECHAMBER ON LIGAND---------"
        frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
        if os.path.exists(self.amberligandfile) and os.path.exists(frcmodfile):
            print "ALREADY HAVE CHARGED MOL2 AND FRCMOD FILES FOR %s" % self.ligand_name
            sys.exit()

        #make new mol2 file for ligand, including charges
        command='%s/bin/antechamber -i %s -fi mol2 -o %s -fo mol2 -c %s -nc %s' % (os.environ['AMBERHOME'], self.ligandfile,  self.amberligandfile, self.charge_method,self.ligand_charge)
        output, err=run_linux_process(command)
        prefix='ante'
        self.check_output(output, err, prefix, type='ante')
        #make frcmod file to catch missing parameters
        prefix='parmchk'
        command='%s/bin/parmchk -i %s -o %s -f mol2' % (os.environ['AMBERHOME'], self.amberligandfile, frcmodfile)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix, type='ante')

        #make amber formattable mol2 file
        #command='%s/bin/antechamber -i %s -fi mol2 -o %s/%s.amber.mol2 -fo mol2 -an y' % (os.environ['AMBERHOME'], self.ligandfile, self.antdir, self.ligand_name)
        #output, err=run_linux_process(command)
        #prefix='mol2convert'
        #self.check_output(output, err, prefix, type='ante')

        # clean up
        #os.system('cd %s' % self.antdir)
        os.system('rm ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF  rm sqm.*')

    def check_output(self, output, err, prefix, type):
        types=['leap', 'ante', 'ptraj', 'md', 'MMGBSA']
        outdirs=[self.leapdir, self.antdir, self.gbdir, self.gbdir, self.gbdir]
        typedict=dict()
        for (t,o) in zip(types, outdirs):
            typedict[type]=o
        dir=typedict[type]
        errors=False
        if 'rror' in err or 'rror' in output:
            errors=True
            if type=='leap' and 'tl_getline()' in err.split('\n')[1]:
                if len(err.split('\n')) < 4:
                    errors=False
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

    def run_leap(self):
        if not os.path.exists(self.leapdir):
            os.mkdir(self.leapdir)
        print "--------------------------------------"
        print "RUNNING LEAP FOR COMPLEX--------------"
        frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
        prefix='cpx'
        reload(amber_file_formatter)
        amber_file_formatter.write_leap(self.leapdir, prefix, self.ligand_name, self.radii, frcmodfile, self.amberligandfile, self.proteinfile, complex=True, gbmin=self.gbmin)
        command='%s/bin/tleap -f %s/%s-%s-leaprc' % (os.environ['AMBERHOME'],
self.leapdir, self.ligand_name, prefix)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix, type='leap')
        # use ante-MMPBSA to make sure topology pieces are all correct
        print "GETTING SEPARATE TOPOLOGIES"
        command='ante-MMPBSA.py -p {0}/{1}-complex.solv.top -c {0}/{1}-complex.solv.crd -s :WAT -r {0}/{1}-protein.top -l {0}/{1}%s-ligand.top
-n :MOL --radii={2}'.format(self.leapdir, self.ligand_name, self.radii)
        return


    def simulation_guts(self, prefix, prmtop, inpcrd):
        # write simulation run input files
        # pass in prefix, prmtop, and inpcrd appropriate for gb vs. explicit
        if 'ligand' in prefix:
            protein_belly=None
        else:
            protein_belly=self.protein_belly
            # store minimized complex name if not just a ligand run
            self.mincpx='%s/%s.rst' % (self.gbdir, prefix)
        print "--------------------------------------"
        if self.gbmin==True:
            print "--------------------------------------"
            print "RUNNING MINIMIZATION WITH GB IMPLICIT-"
            amber_file_formatter.write_simulation_input(md=False, dir=self.gbdir, prefix=prefix, gbmin=self.gbmin, gb_model=self.gb_model, protein_belly=protein_belly, maxcycles=self.maxcycles)
        else:
            print "RUNNING MINIMIZATION WITH EXPLICIT----"
            amber_file_formatter.write_simulation_input(md=False, dir=self.gbdir, prefix=prefix, protein_belly=protein_belly, maxcycles=self.maxcycles)

        command=self.get_simulation_commands(prefix, prmtop, inpcrd)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix=prefix, type='md')
        # if ligand minimization, totally skip MD
        if 'ligand' in prefix:
            pass
        elif self.md==True:
            print "--------------------------------------"
            inpcrd=self.mincpx
            if self.gbmin==True:
                print "RUNNING MD SIMULATION WITH GB IMPLICIT"
                amber_file_formatter.write_simulation_input(md=True, dir=self.gbdir, prefix=prefix,  gbmin=self.gbmin, gb_model=self.gb_model, protein_belly=protein_belly)
            else:
                print "RUNNING MD SIMULATION WITH EXPLICIT---"
                amber_file_formatter.write_simulation_input(md=True, dir=self.gbdir, prefix=prefix,  protein_belly=protein_belly)
            command=self.get_simulation_commands(prefix, prmtop, inpcrd)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='md', type='md')
        else:
            pass
        return
        

    def run_cpx_simulation(self):
        # set intput variables and commands
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
        if self.ligand_restraints==True and self.protein_radius!=None:
            protein_belly=get_belly_restraints(self.protein_radius, prmtop, inpcrd, ligand_restraints=True)
        elif self.ligand_restraints!=True and self.protein_radius!=None:
            protein_belly=get_belly_restraints(self.protein_radius, prmtop, inpcrd, ligand_restraints=False)
        else:
            protein_belly=None
        self.protein_belly=protein_belly
        simulation_guts(self, prefix, prmtop, inpcrd)
        if self.md==True:
            simulation_guts(self, mdprefix, prmtop, inpcrd)
        return

    def run_ligand_strain(self):
        # extract ligand from minimized complex structure with cpptraj
        print "--------------------------------------"
        print "SETTING UP LIGAND STRAIN CALC---------"
        filename='%s/getligand.ptraj' % self.gbdir
        if self.gbmin==True:
            self.minligandcpx='%s/%s-ligand-cpxmin.rst' % (self.leapdir, self.ligand_name)
            prmtop='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
            # now strip with ptraj to get an restart file for MMPBSA
            amber_file_formatter.write_ptraj_strip(filename, self.mincpx, self.minligandcpx)
            command='cpptraj {0}/{1}-complex.top {2}'.format(self.leapdir, self.ligand_name, filename)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='strip-ligand', type='ptraj')
            prefix='mingb-ligand'
        else:
            # rebuild ligand topology only if using explicit solvent
            # use minimized ligand outconf mol2 as input structure
            outconf='%s/%s-ligand-cpxmin.solv.mol2' % (self.leapdir, self.ligand_name)
            frcmodfile='%s/%s.frcmod' % (self.antdir, self.ligand_name)
            amber_file_formatter.write_leap(dir=self.leapdir, prefix=leap_prefix, ligand_name=self.ligand_name, radii=self.radii, frcmodfile=frcmodfile, newligandfile=outconf, complex=False, gbmin=False)
            command='{0}/bin/tleap -f {1}/{2}-{3}-leaprc'.format(os.environ['AMBERHOME'], self.leapdir, self.ligand_name, prefix)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='ligsolv', type='leap')
            print "RAN LEAP for RESOLVATING LIGAND"
            prmtop='%s/%s-ligand.solv.top' % (self.leapdir, self.ligand_name)
            # now strip with ptraj
            amber_file_formatter.write_ptraj_strip(filename, self.mincpx, outconf)
            command='cpptraj {0}/{1}-complex.solv.top {2}'.format(self.leapdir, self.ligand_name, filename)
            output, err=run_linux_process(command)
            self.check_output(output, err, prefix='strip-ligand', type='ptraj')
            prefix='min-ligand'
        # minimize ligand in solution (or with GB solvent)
        # see diff in complex energy and solvated energy
        simulation_guts(prefix, prmtop, outconf)
        print "MINIMIZED LIGAND"
        self.run_mmgbsa(complex=False)
        print "MMGB CALC FINISHED ON LIGAND"
        return

    def mmgbsa_guts(self, prefix, start, finish, solvcomplex, complex, traj, interval=1, protein=None, ligand=None):
        inputfile='%s/%s-mmgb.in' % (self.gbdir, prefix)
        amber_file_formatter.write_mmgbsa_input(inputfile, self.gb_model, start, interval, finish)
        # use MMGBSA.py in Amber14 to run MMGB free energy difference calcs for complex
        if protein!=None and ligand!=None:
            print "--------------------------------------"
            print "RUNNING COMPLEX MMGBSA CALC-----------"
            command='{0}/bin/MMPBSA.py -i {1} -o {2}/{3}-{4}-FINAL_MMPBSA.dat -sp {5} -cp {6} -rp {7} -lp {8} -y {9}'.format(os.environ['AMBERHOME'], inputfile, self.gbdir, self.ligand_name, prefix, solvcomplex, complex, protein, ligand, traj)
        else:
            print "--------------------------------------"
            print "RUNNING LIGAND MMGBSA CALC------------"
            command='{0}/bin/MMPBSA.py -i {1} -o {2}/{3}-{4}-FINAL_MMPBSA.dat -sp {5} -cp {6} -y {7}'.format(os.environ['AMBERHOME'], inputfile, self.gbdir, self.ligand_name, prefix, solvcomplex, complex, traj)
        output, err=run_linux_process(command)
        self.check_output(output, err, prefix, type='MMGBSA')
        return

    def run_mmgbsa(self, complex=True):
        #setup files for MMGBSA.py in Amber14 to run MMGB free energy difference calcs for complex
        inputfile='%s/mmgb.in' % self.gbdir
        if complex==True:
            prefix='cpx'
            start=0
            interval=1
            if self.md==True:
                finish=-1
                if self.gbmin==True:
                    traj='%s/gbmd.mdcrd' % (self.gbdir)
                    solvcomplex=None
                else:
                    traj='%s/md.mdcrd' % (self.gbdir)
                    solvcomplex='%s/%s-complex.solv.top' % (self.leapdir, self.ligand_name)
            else:
                finish=1
                if self.gbmin==True:
                    traj='%s/gbmin.rst' % (self.gbdir)
                    solvcomplex=None
                else:
                    traj='%s/min.rst' % (self.gbdir)
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
                solvcomplex='%s/%s-ligand.top' % (self.leapdir, self.ligand_name)
                initial_traj=self.minligandcpx
                final_traj='%s/gbmin-ligand.rst' % (self.gbdir)
            else:
                solvcomplex='%s/%s-ligand.solv.top' % (self.leapdir, self.ligand_name)
                initial_traj='%s/%s-ligand.solv.crd' % (self.leapdir, self.ligand_name)
                final_traj='%s/min-ligand.rst' % (self.gbdir)
            # first get initial GB energy of ligand in complex
            prefix='ligcpx'
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, initial_traj, interval=1)
            # next GB energy of ligand minimized in solution
            prefix='ligsolv'
            self.mmgbsa_guts(prefix, start, finish, solvcomplex, complex, final_traj, interval=1)
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
                    if 'EEL' in line:
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
        keyorder=['MMGB-strain', 'MMGB', 'strain', 'vdW', 'eel_inter', 'eel/EGB', 'EGB', 'E_surf', 'E_lig'] 
        entry=''.join(['%s\t' % x for x in keyorder])
        entry=''.join(['name\t', entry, '\n'])
        ohandle.write(entry)
        print entry
        for ligand in sorted_ligands:
            all_values[ligand]['MMGB-strain']+=all_values[ligand]['MMGB']+all_values[ligand]['strain']
            name='%s\t' % ligand
            entry=''.join(['%0.3f\t' % round(float(all_values[ligand][x]), 3) for x in keyorder])
            entry=''.join([name, entry, '\n'])
            ohandle.write(entry)
            print entry


