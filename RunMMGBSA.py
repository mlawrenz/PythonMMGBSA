import PythonMMGBSA
import time
import optparse


def main(proteinfile, ligandfile, ligand_charge, gb_model, maxcycles, md, steps=0, gbmin=False, gpu=False):

    mol=PythonMMGBSA.ambermol(proteinfile, ligandfile, gb_model=gb_model,
ligand_charge=ligand_charge, gpu=gpu, protein_radius=0.1, maxcycles=maxcycles, md=md, md_steps=steps, gbmin=gbmin)
    start=time.time()
    mol.run_antechamber()
    end=time.time()
    print "antechamber elapased time %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_leap()
    end=time.time()
    print "leap elapsed time %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_cpx_simulation() 
    end=time.time()
    print "cpx sim elapsed time %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_ligand_strain()
    end=time.time()
    print "lig strain elapsed time %s sec" % str(round(end-start, 2))
    start=time.time()
    mol.run_mmgbsa(complex=True)
    end=time.time()
    print "cpx mmgbsa elapsed time %s sec" % str(round(end-start, 2))
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
    parser.add_option('-s','--steps',dest='steps',type='string', help='MD simulation steps (2 fs)')
    parser.add_option('-u', action="store_true", dest="gpu", help="using flag -p will run a GPU MD simulation")
    parser.add_option('-m', action="store_true", dest="md", help="using flag -m will run a MD simulation")
    parser.add_option('-g', action="store_true", dest="gbmin", help="using flag -g will run implicit GB solvent instead of explicit solvent simulations")


    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    if options.maxcycles==None:
        options.maxcycles=50000
    if options.steps==None:
        options.steps=100000
    if options.md==True:
        if options.gbmin==True:
            if options.gpu==True:
                main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=True, steps=options.steps, gpu=True, gbmin=True)
            else:
                main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=True, steps=options.steps, gpu=False, gbmin=True)
        else:
            if options.gpu==True:
                main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=True, steps=options.steps, gpu=True)
            else:
                main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=True, steps=options.steps, gpu=False)
    else:
        if options.gbmin==True:
            main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=False, gbmin=True)
        else:
            main(proteinfile=options.proteinfile, ligandfile=options.ligandfile, ligand_charge=options.ligand_charge, gb_model=options.gb_model, maxcycles=options.maxcycles, md=False)

