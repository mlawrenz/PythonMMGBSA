import os, time
import PythonMMGBSA

# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(epilog='''This program will calculate binding free
                        energies using end-state free energy methods.''', formatter_class=ArgumentDefaultsHelpFormatter)

# for turning on timer
parser.add_argument('-time', action="store_true", dest="time", help="using -time will turn on timing of actions")

#always need these user inputs
group = parser.add_argument_group('Necessary user input',)
group.add_argument('-prot','--protfile',dest='protfile',  help='protein PDB file')
group.add_argument('-mol2','--ligfile',dest='ligfile',  help='ligand MOL2 file, named MOL')
group.add_argument('-netc','--netcharge',dest='ligcharge',  help='total net charge on ligand')
# opions with defaults
group = parser.add_argument_group('Options with defaults',)
group.add_argument('-gb','--gbmodel',dest='gbmodel',  help='MMGB model version in AMBER', default=8)
group.add_argument('-prad','--proteinrad',dest='prot_radius',  help='distance around ligand allowed to move in protein', default=0.1)
group.add_argument('-ligr', action="store_true", dest="ligrestraint", help="using flag -ligr will restrain ligand atoms by k=5.0")
group.add_argument('-drms','--drms',dest='drms',  help='max rmsd of energy gradient ', default=0.1)
group.add_argument('-maxcyc','--maxcycles',dest='maxcycles',  help='min maxcycles', default=50000)
group.add_argument('-im', action="store_true", dest="gbmin", help="using flag -im will run implicit GB solvent instead of explicit solvent simulations")
group.add_argument('-gpu', action="store_true", dest="gpu", help="using flag -p will run a GPU MD simulation", default=False)
# options specific for md
group = parser.add_argument_group('Options that turn on MD')
group.add_argument('-md', action="store_true", dest="md", help="using flag -m will run a MD simulation")
group.add_argument('-mdsteps','--mdsteps',dest='mdsteps',  help='MD simulation steps (2 fs)')


def main(args):
    mol=PythonMMGBSA.ambermol(protfile=args.protfile, ligfile=args.ligfile, ligcharge=args.ligcharge, gbmodel=args.gbmodel, prot_radius=args.prot_radius, ligrestraint=args.ligrestraint, maxcycles=args.maxcycles, drms=args.drms, gbmin=args.gbmin, md=args.md, mdsteps=args.mdsteps)
    #mol.run_antechamber()
    #mol.run_leap()
    #mol.run_cpx_simulation() 
    #mol.run_ligand_strain()
    #mol.run_mmgbsa(complex=True)
    mol.print_table()


def time_main(args):
    print "timer turned on!"
    mol=PythonMMGBSA.ambermol(protfile=args.protfile, ligfile=args.ligfile, ligcharge=args.ligcharge, gbmodel=args.gbmodel, prot_radius=args.prot_radius, ligrestraint=args.ligrestraint, maxcycles=args.maxcycles, drms=args.drms, gbmin=args.gbmin, md=args.md, mdsteps=args.mdsteps)
    start=time.time()
    mol.run_antechamber()
    end=time.time()
    print "antechamber ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_leap()
    end=time.time()
    print "leap ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_cpx_simulation() 
    end=time.time()
    print "cpx simulation ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_ligand_strain()
    end=time.time()
    print "ligand strain calc ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_mmgbsa(complex=True)
    end=time.time()
    print "cpx mmgbsa ran in %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.print_table()


if __name__=="__main__":	
    args = parser.parse_args()
    if args.time==True:
        time_main(args)
    else:
        main(args)

