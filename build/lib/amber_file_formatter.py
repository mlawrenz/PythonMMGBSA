import sys
import subprocess
import pickle
import optparse
import pylab
import numpy
import glob
import os
from os.path import split

def write_leap(dir, prefix, ligand_name, radii, frcmodfile, newligfile,
protfile=None, complex=True, implicit=False):
        fhandle=open('{0}/{1}-{2}-leaprc'.format(dir, ligand_name, prefix), 'w')
        # write headers common to all
        fhandle.write('''\
source leaprc.ff14SB
source leaprc.gaff
loadAmberParams frcmod.ionsjc_tip3p

loadAmberParams %s
''' % frcmodfile)
        if complex==True:
            # write complex info 
            fhandle.write('''\
mol=loadmol2 {0}
prot=loadpdb {1}
complex = combine {{prot mol}}

set default PBradii {2}
saveAmberParm mol {4}/{3}-ligand.top {4}/{3}-ligand.crd
saveAmberParm prot {4}/{3}-protein.top {4}/{3}-protein.crd


saveAmberParm complex {4}/{3}-complex.top {4}/{3}-complex.crd
savepdb complex {4}/{3}-complex.pdb
'''.format(newligfile,protfile, radii, ligand_name, dir))
            if implicit==False:
                fhandle.write('''\
solvateOct complex TIP3PBOX 12.0
addIons complex Na+ 0
saveAmberParm complex {1}/{0}-complex.solv.top {1}/{0}-complex.solv.crd
savepdb complex {1}/{0}-complex.solv.pdb
quit'''.format(ligand_name, dir))
            else:
                pass # don't solvate system
        else:
            # only need to rewrite ligand topo if solvating it
            fhandle.write('''\
mol=loadmol2 {0}
solvateOct mol TIP3PBOX 12.0
set default PBradii {3}
saveAmberParm mol {1}/{2}-ligand.solv.top {1}/{2}-ligand.solv.crd
savepdb mol {1}/{2}-ligand.solv.amber.pdb'''.format( newligfile, dir, ligand_name, radii))
        fhandle.close()

# previously used belly dynamics, but igb does not work with belly
def add_belly(fhandle, protein_belly):
    fhandle.write('''\
  ibelly=1, bellymask = ":%s"''' % protein_belly)
    return fhandle

def add_restraints(fhandle, restraint_atoms, restraint_k=10.0):
    fhandle.write('''\
  ntr = 1, 
  restraintmask = ":{0}", restraint_wt = {1},'''.format(restraint_atoms, restraint_k))
    return fhandle

def write_simulation_input(md, dir, prefix, gbmodel=0, restraint_atoms=None, \
restraint_k=10.0, maxcycles=50000, drms=0.1, steps=100000, mdseed=-1):
    heatsteps=25000
    totalsteps=int(heatsteps)+int(steps)
    filename='%s/%s.in' % (dir, prefix)
    if gbmodel!=0:
        periodic=0 # no box if gb is used
    else:
        periodic=1
    fhandle=open(filename, 'w')
    if md==False:
        fhandle.write('''\
minimization
  &cntrl
  imin = 1, maxcyc = {0}, ntmin = 1,
  ncyc   = 100, 
  ntx = 1, ntc = 1, ntf = 1,
  ntb = {1}, ntp = 0,
  ntwx = 1000, ntwe = 0, ntpr = 1000,
  igb={2},
  drms   = {3},
  cut = 10.0,'''.format(maxcycles,periodic, gbmodel, drms))
        if restraint_atoms!=None:
            fhandle=add_restraints(fhandle, restraint_atoms, restraint_k)
        fhandle.write('&end\n')
    else:
        fhandle.write('''\
nvt equilibration with Langevin therm, SHAKE Hbonds
  &cntrl
  imin = 0, ntx = 1, irest = 0, nstlim = {0},
  ntt=3,  gamma_ln=1.0, temp0 = 298.15, tempi = 0, ig = {1},
  ntc = 2, ntf = 2, dt = 0.002,
  ntb = {2}, ntp = 0,
  ntwx = 1000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 10.0, iwrap = 1,
  nscm = 100,
  igb={3},
  nscm = 100,'''.format(totalsteps, mdseed, periodic, gbmodel))
        if restraint_atoms!=None:
            fhandle=add_restraints(fhandle, restraint_atoms, restraint_k)
        fhandle.write('&end\n')
        fhandle.write('&wt\n')
        fhandle.write('''
step1=0, istep2={0}, value1=0, value2=300,
&end
&wt 
type='TEMP0', istep1={0}, istep2={1}, value1=300, value2=300,
&end'''.format(heatsteps, totalsteps))
        fhandle.close()
    return

def write_mmgbsa_input(filename, model, start, interval, finish=100000000):
    # script will redice final frame to total frames if really high
    fhandle=open(filename, 'w')
    fhandle.write('''\
run GB 
&general
 startframe={0}, interval={1}, endframe={2},
 keep_files=0,
 /
&gb
 igb={3}, saltcon=0.15,
 /'''.format( start, interval, finish, model))
    fhandle.close()

def write_ptraj_strip(filename, inconf, outconf):
    filetypes=['mol2', 'rst', 'crd', 'pdb']
    mapper=dict()
    mapper['mol2']='mol2'
    mapper['rst']='restart'
    mapper['crd']='inpcrd'
    mapper['pdb']='pdb'
    files=dict()
    types=dict()
    files['in']=inconf
    files['out']=outconf
    for option in files.keys():
        base=os.path.basename(files[option])
        suffix=base.split('.')[-1]
        types[option]=mapper[suffix]
        if suffix not in filetypes:
            print "FILETYPE NOT SUPPORTED FOR PTRAJ: %s" % types[option]
            sys.exit()
    fhandle=open(filename, 'w')
    fhandle.write('''\
trajin {0} {1}
strip !:MOL
trajout {2} {3}'''.format(inconf, types['in'],  outconf, types['out'])) 
    fhandle.close()

