import PythonMMGBSA
import optparse



def main(proteinfile, ligandfile):
    mol=PythonMMGBSA.ambermol(proteinfile, ligandfile)
    mol.run_antechamber()
    mol.run_leap()
    mol.run_simulation()
    mol.run_mmgbsa()


def parse_cmdln():
    import os
    parser=optparse.OptionParser()
    parser.add_option('-p','--proteinfile',dest='proteinfile',type='string')
    parser.add_option('-l','--ligandfile',dest='ligandfile',type='string')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__=="__main__":	
    (options,args)=parse_cmdln()
    main(proteinfile=options.proteinfile, ligandfile=options.ligandfile)
