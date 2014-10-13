import PythonMMGBSA
import optparse


def main(proteinfile, ligandfile, ligand_charge, gb_model, md, gbmin=False, maxcycles=50000):

    mol=PythonMMGBSA.ambermol(proteinfile, ligandfile, gb_model=gb_model, ligand_charge=ligand_charge, gpu=False, protein_radius=0.1, maxcycles=maxcycles, md=md, gbmin=gbmin)
    mol.run_antechamber()
    mol.run_leap()
    mol.run_cpx_simulation() 
    mol.run_ligand_strain()
    mol.run_mmgbsa(complex=True)
    # run after generating all results to get a table
    mol.print_table()


def parse_cmdln():
    import os
    parser=optparse.OptionParser()
    parser.add_option('-p','--proteinfile',dest='proteinfile',type='string', help='protein PDB file')
    parser.add_option('-l','--ligandfile',dest='ligandfile',type='string', help='ligand MOL2 file, named MOL')
    parser.add_option('-c','--charge',dest='ligand_charge',type='string', help='total charge on ligand')
    parser.add_option('-b','--gb_model',dest='gb_model',type='string', help='MMGB model version in AMBER')
    parser.add_option('-y','--maxcycles',dest='maxcycles',type='string', help='min maxcycles')
    parser.add_option('-m', action="store_true", dest="md", help="using flag -m will run a MD simulation")
    parser.add_option('-g', action="store_true", dest="gbmin", help="using flag -g will run implicit GB solvent instead of explicit solvent simulations")


    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    if options.maxcycles==None:
        options.maxcycles=50000
    if options.md==True:
        if options.gbmin==True:
            main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=True, gbmin=True)
        else:
            main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=True)

    else:
        if options.gbmin==True:
            main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=False, gbmin=True)
        else:
            main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=False)



