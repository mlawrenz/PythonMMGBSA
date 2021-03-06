import os, time
import shutil
import sys
import amber_file_formatter
import PythonMMGBSA

# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(epilog='''This program will calculate binding free
                        energies using end-state free energy methods.''', formatter_class=ArgumentDefaultsHelpFormatter)

#always need these user inputs
group = parser.add_argument_group('Necessary user input',)
group.add_argument('-mol2','--ligfile',dest='ligfile',  help='ligand MOL2 file, named MOL')
group.add_argument('-netc','--netcharge',dest='ligcharge',  help='total net charge on ligand')
group.add_argument('-im', action="store_true", dest="gbmin", help="using flag -im will run implicit GB solvent instead of explicit solvent simulations", default=False)
group.add_argument('-gb','--gbmodel',dest='gbmodel',  help='MMGB model version in AMBER', default=1)


def antechamber(antdir, ligfile, ligcharge):
    origdir=os.getcwd()
    os.chdir(antdir)
    ligandname=os.path.basename(ligfile).split('.mol2')[0]
    amberligfile='%s.amber.mol2' % ligandname
    charge_method='bcc'
    command='%s/bin/antechamber -i %s/%s -fi mol2 -o %s.amber.pdb -fo pdb' % (os.environ['AMBERHOME'], origdir, ligfile, ligandname)
    output, err=PythonMMGBSA.run_linux_process(command)
    if 'rror' in err or 'rror' in output:
        print "ERRORS"
        print err, output
    command='{0}/bin/antechamber -i {1}.amber.pdb -fi pdb -o {1}.prep -fo prepi -c {2} -nc {3}'.format(os.environ['AMBERHOME'], ligandname,charge_method,ligcharge) 
    print command
    output, err=PythonMMGBSA.run_linux_process(command)
    if 'rror' in err or 'rror' in output:
        print "ERRORS"
        print err, output
    command='{0}/bin/parmchk -i {1}.prep -o {1}.frcmod -f prepi'.format(os.environ['AMBERHOME'], ligandname)
    output, err=PythonMMGBSA.run_linux_process(command)
    if 'rror' in err or 'rror' in output:
        print "ERRORS"
        print err, output
        sys.exit()
    os.system('rm ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF  rm sqm.*')
    print "RAN ANTECHAMBER"
    os.chdir(origdir)
    return

def writeleap(tmpdir, antdir, ligandname, radii, gbmin):
    fhandle=open('%s/leaprc' % tmpdir, 'w')
    fhandle.write('''\
source leaprc.ff14SB
source leaprc.gaff
loadAmberParams frcmod.ionsjc_tip3p

loadAmberParams {0}/{1}.frcmod
loadamberprep {0}/{1}.prep

mol=loadpdb {3}/{1}.amber.pdb
set default PBradii {2}
saveAmberParm mol {3}/{1}.top {3}/{1}.crd
'''.format(antdir, ligandname, radii, tmpdir))
    if gbmin==False:
        fhandle.write('''\
solvateOct mol TIP3PBOX 14.0
saveAmberParm mol {0}/{1}.solv.top {0}/{1}.solv.crd'''.format(tmpdir, ligandname))
    fhandle.write('\nquit')
    return

def main(args):
    ligfile=args.ligfile
    ligcharge=args.ligcharge
    gbmin=args.gbmin
    gbmodel=args.gbmodel
    radii=PythonMMGBSA.get_pbbond_radii(gbmodel)
    charge_method='bcc'
    origdir=os.getcwd()
    fullligandname=os.path.basename(ligfile).split('.mol2')[0]
    ligandname=os.path.basename(ligfile).split('-')[0]
    tmpdir='%s-tmp' % ligandname
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    # RUN ANTECHAMBER IF YOU NEED TO
    antdir='{0}-antechamber-output/'.format(ligandname)
    if not os.path.exists(antdir):
        os.mkdir(antdir)
    if os.path.exists('%s/%s.prep' % (antdir, ligandname)):
        print "already have antechamber files for %s" % ligandname
        # get PDB though
        command='%s/bin/antechamber -i %s -fi mol2 -o %s/%s.amber.pdb -fo pdb' % (os.environ['AMBERHOME'], ligfile, tmpdir, ligandname)
        output, err=PythonMMGBSA.run_linux_process(command)
        if 'rror' in err or 'rror' in output:
            print "ERRORS"
            print err, output
            sys.exit()
    else:
        antechamber(antdir, ligfile, ligcharge)
    if gbmin==True:
        prmtop='%s/%s.top' % (tmpdir, ligandname)
        inpcrd='%s/%s.crd' % (tmpdir, ligandname)
    else:
        prmtop='%s/%s.solv.top' % (tmpdir, ligandname)
        inpcrd='%s/%s.solv.crd' % (tmpdir, ligandname)
    writeleap(tmpdir, antdir, ligandname, radii, gbmin)
    command='{0}/bin/tleap -f {1}/leaprc'.format(os.environ['AMBERHOME'], tmpdir)
    output, err=PythonMMGBSA.run_linux_process(command)
    if 'rror' in err or 'rror' in output:
        print "ERRORS"
        print err, output
        sys.exit()
    print "RAN LEAP"
    # DO MINIMIZATION
    if gbmin==True:
        prefix='gbmin'
    else:
        prefix='min'
    amber_file_formatter.write_simulation_input(md=False, dir=tmpdir, prefix=prefix,  gbmin=gbmin, gbmodel=gbmodel, drms=0.001)
    command=PythonMMGBSA.get_simulation_commands(prefix, prmtop, inpcrd, tmpdir, gpu=False, restrain=False, nproc=16, mdrun=False)
    output, err=PythonMMGBSA.run_linux_process(command)
    if 'rror' in err or 'rror' in output:
        print "ERRORS"
        print err, output
        sys.exit()
    print "done with sim"
    amber_file_formatter.write_mmgbsa_input('%s/mmgb.in' % tmpdir, gbmodel, start=0, interval=1, finish=1)
    command='MMPBSA.py -i {0}/mmgb.in -o {0}/ENERGY.dat -sp {0}/{2}.solv.top -cp {0}/{2}.top -y {0}/min.rst'.format(tmpdir, fullligandname, ligandname)
    output, err=PythonMMGBSA.run_linux_process(command)
    if 'rror' in err or 'rror' in output:
        print "ERRORS"
        print err, output
        sys.exit()
    print "ran MMPBSA on conf"
    command='echo {0} >> {1}-ENERGIES.OUT'.format(fullligandname, ligandname)
    output, err=PythonMMGBSA.run_linux_process(command)
    command='sed "0,/GENER/d" < {0}/ENERGY.dat >> {1}-ENERGIES.OUT'.format(tmpdir, ligandname)
    output, err=PythonMMGBSA.run_linux_process(command)
    shutil.rmtree(tmpdir)


if __name__=="__main__":	
    args = parser.parse_args()
    main(args)

