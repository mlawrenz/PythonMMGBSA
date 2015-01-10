import os, time
import PythonMMGBSA

# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(epilog='''This program will calculate binding free
                        energies using end-state free energy methods.''', formatter_class=ArgumentDefaultsHelpFormatter)

# true/false options
parser.add_argument('-time', action="store_true", dest="time", help="Using -time will turn on timing of actions.")
parser.add_argument('-keepfiles', action="store_true", dest="keepfiles", help="Using -keepfiles will keep all intermediate files.")

# always need these user inputs
group = parser.add_argument_group('Necessary user input:',)
group.add_argument('-jname','--jobname',dest='jobname',  help='Label prefix for output directory.')
group.add_argument('-prot','--protfile',dest='protfile',  help='Protein PDB file.')
group.add_argument('-mol2','--ligfile',dest='ligfile',  help='Ligand MOL2 file.\
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

# options specific for md
group = parser.add_argument_group('Options that pertain to MD:')
group.add_argument('-md', action="store_true", dest="md", help="Using flag -m will run a MD simulation after minimization.")
group.add_argument('-nproc','--nproc',dest='nproc',  help='Number of processors to run MPI processes.', default=8)
group.add_argument('-mdsteps','--mdsteps',dest='mdsteps',  help='MD simulation time steps (1 step=2 fs).', default=100000)
group.add_argument('-gpu', action="store_true", dest="gpu", help="Using flag -p will run a GPU MD simulation.", default=False)


def main(args):
    mol=PythonMMGBSA.ambermol(jobname=args.jobname, protfile=args.protfile, ligfile=args.ligfile, ligcharge=args.ligcharge, gbmodel=args.gbmodel, prot_radius=args.prot_radius, maxcycles=args.maxcycles, drms=args.drms, implicit=args.implicit, md=args.md, mdsteps=args.mdsteps, nproc=args.nproc, gpu=args.gpu)
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
    mol=PythonMMGBSA.ambermol(jobname=args.jobname, protfile=args.protfile, ligfile=args.ligfile, ligcharge=args.ligcharge, gbmodel=args.gbmodel, prot_radius=args.prot_radius, maxcycles=args.maxcycles, drms=args.drms, implicit=args.implicit, md=args.md, mdsteps=args.mdsteps, nproc=args.nproc, gpu=args.gpu)
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
    if args.time==True:
        time_main(args)
    else:
        main(args)

