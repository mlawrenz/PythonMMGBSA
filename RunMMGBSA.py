import os, time
import sys
import PythonMMGBSA

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter

parser = ArgumentParser(description='''
Written by Morgan Lawrenz, 2014
This program will calculate binding free energies for protein-ligand complexes using
the end-state MMGBSA free energy method, and includes a ligand strain penalty.

The program calls the python implementation of AMBER MMPBSA described in this paper:
Miller III, B. R., McGee Jr., T. D., Swails, J. M. Homeyer, N. Gohlke, H. and 
Roitberg, A. E. J. Chem. Theory Comput., 2012, 8 (9) pp 3314--3321 

Source code and README file in /common/compchem/mlawrenz/PythonMMGBSA/.

This script runs through 6 main modules for the workflow:

    mol=PythonMMGBSA.ambermol(args) ----> mol object is initialized with input
    mol.run_antechamber() ----> generates ligand AMBER parameters
    mol.run_leap() ----> builds complete (solvated) complex topology file
    mol.run_cpx_simulation() ----> runs minimization or MD simulation
    mol.run_ligand_strain() ---->  computes ligand strain energy
    mol.run_mmgbsa(complex=True) ---->  computes MMGBSA free energy of binding
    mol.print_table() ---->  prints results.tbl in the output directory
    mol.clean() ----> removes all but final MMGB output and relevant structure files

You can customize this script by only calling select steps.''', epilog='''
**The script does not automatically launch LSF jobs**
See /common/compchem/mlawrenz/PythonMMGBSA/example-submission-script.sh or
bsub -J mmgb -n nproc -q queue -R span[hosts=1] RunMMGBSA.py [ args ]''', formatter_class=RawTextHelpFormatter)




# true/false options
parser.add_argument('-time', action="store_true", dest="time", help="Using -time will turn on timing of actions.")
parser.add_argument('-keepfiles', action="store_true", dest="keepfiles", help="Using -keepfiles will keep all intermediate files.")

# always need these user inputs
group = parser.add_argument_group('Necessary user input:',)
group.add_argument('-jname','--jobname',dest='jobname',  help='Label prefix for output directory.')
group.add_argument('-prot','--protfile',dest='protfile',  help='Protein PDB file.')
group.add_argument('-mol2','--ligfile',dest='ligfile',  help='Ligand MOL2 file. \
IMPORTANT: name ligand residue MOL and ensure ligand coordinates are the bound pose in the pocket of the separately provided protfile')
group.add_argument('-netc','--netcharge',dest='ligcharge',  help='Total net charge on ligand.')


# options with defaults
group = parser.add_argument_group('Options with defaults:',)
group.add_argument('-gb','--gbmodel',dest='gbmodel',  help='MMGB model version in AMBER.', default=1)
group.add_argument('-prad','--proteinrad',dest='prot_radius',  help='Distance \
around ligand to select moveable protein residues. Can set to None in main() for no restraints.', default=0.1)
group.add_argument('-drms','--drms',dest='drms',  help='Max rmsd of energy gradient during minimization.', default=0.1)
group.add_argument('-maxcyc','--maxcycles',dest='maxcycles',  help='Max cycles of minimization.', default=50000)
group.add_argument('-im', action="store_true", dest="implicit", help="Using flag \
-im will run implicit GB solvent: NOT necessarily faster for MD. Please benchmark before using implicit instead of explicit for MD.")
group.add_argument('-mdseed','--mdseed',dest='mdseed',  help='MD seed for random number generator; default ensures unique seed per simulation.', default=-1)

# options specific for md
group = parser.add_argument_group('Options that pertain to MD:')
group.add_argument('-md', action="store_true", dest="md", help="Using flag -m will run a MD simulation after minimization.")
group.add_argument('-nproc','--nproc',dest='nproc',  help='Number of processors to run MPI processes.', default=8)
group.add_argument('-mdsteps','--mdsteps',dest='mdsteps',  help='MD simulation time steps (1 step=2 fs).', default=100000)
group.add_argument('-gpu', action="store_true", dest="gpu", help="Using flag -gpu will run a GPU MD simulation.", default=False)


def check_environment():
    # check for mpirun, amber14
    if 'amber14' not in os.environ['AMBERHOME']:
        print "NEED TO RESET AMBERHOME to /common/compchem/src/amber14/"
        sys.exit()
    path_items=os.environ['PATH'].split(':')
    amber_paths=[]
    mpi=False
    amber=False
    if '/common/compchem/src/amber14/lib' not in os.environ['LD_LIBRARY_PATH']:
        print "NEED AMBER14 LIB IN LD_LIBRARY_PATH: /common/compchem/src/amber14/lib"
        sys.exit()
    for path in path_items:
        if 'amber14-mpirun-build' in path:
            print "HAVE WORKING MPIRUN: ", path
            mpi=True
        if 'mpich-3.1.2-build' in path:
            print "HAVE WORKING MPIRUN: ", path
            mpi=True
        if 'amber' in path:
            amber=True
            amber_paths.append(path)
    
    if amber==False:
        print "AMBER14 NEEDS TO BE IN PATH: /common/compchem/src/amber14/bin"
        sys.exit()
    if mpi==False:
        print "NEED WORKING MPIRUN in path: /common/compchem/src/amber14-mpirun-build/"
        sys.exit()
    if len(amber_paths) > 1:
        if 'amber14' not in amber_paths[0]:
            print "TWO AMBER DISTRIBUTIONS IN PATH"
            print "AMBER14 NEEDS TO BE FIRST: /common/compchem/src/amber14/bin"
            sys.exit()
        else:
            print "AMBER14 FOUND"
    else:
        if 'amber14' not in amber_paths[0]:
            print "AMBER14 NEEDS TO BE IN PATH: /common/compchem/src/amber14/bin"
        else:
            print "AMBER14 FOUND"
    return 
    

def main(args):
    mol=PythonMMGBSA.ambermol(jobname=args.jobname, protfile=args.protfile, \
ligfile=args.ligfile, ligcharge=args.ligcharge, gbmodel=args.gbmodel, \
prot_radius=args.prot_radius, maxcycles=args.maxcycles, drms=args.drms, \
implicit=args.implicit, md=args.md, mdsteps=args.mdsteps, mdseed=args.mdseed, nproc=args.nproc, gpu=args.gpu)
    mol.run_antechamber()
    mol.run_leap()
    mol.run_cpx_simulation() 
    mol.run_ligand_strain()
    mol.run_mmgbsa(complex=True)
    mol.print_table()
    if args.keepfiles!=True:
        mol.clean()


def time_main(args):
    print "Timer turned on!"
    mol=PythonMMGBSA.ambermol(jobname=args.jobname, protfile=args.protfile, \
ligfile=args.ligfile, ligcharge=args.ligcharge, gbmodel=args.gbmodel,prot_radius=args.prot_radius,\
 maxcycles=args.maxcycles, drms=args.drms, implicit=args.implicit, md=args.md, mdsteps=args.mdsteps, \
mdseed=args.mdseed, nproc=args.nproc, gpu=args.gpu)
    start=time.time()
    mol.run_antechamber()
    end=time.time()
    print "Antechamber ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_leap()
    end=time.time()
    print "Leap ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_cpx_simulation() 
    end=time.time()
    print "Cpx simulation ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_ligand_strain()
    end=time.time()
    print "Ligand strain calc ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_mmgbsa(complex=True)
    end=time.time()
    print "Cpx mmgbsa ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.print_table()
    if args.keepfiles!=True:
        mol.clean()


if __name__=="__main__":	
    args = parser.parse_args()
    check_environment()
    if args.time==True:
        time_main(args)
    else:
        main(args)

