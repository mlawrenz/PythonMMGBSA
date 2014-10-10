import sys
import subprocess
import pickle
import optparse
import pylab
import numpy
import glob
import os
from os.path import split

def write_leap(dir, prefix, ligand_name, radii, frcmodfile, newligandfile, proteinfile=None, complex=True, gbmin=False):
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
'''.format((newligandfile,proteinfile, radii, ligand_name, dir))
            if gbmin==False:
                fhandle.write('''\
solvateOct complex TIP3PBOX 14.0
#addIons complex Na+ 0
saveAmberParm complex {1}/{0}-complex.solv.top {1}/{0}-complex.solv.crd
quit'''.format(ligand_name, dir))
            else:
                pass # don't solvate system
        else:
            # only need to rewrite ligand topo if solvating it
            fhandle.write('''\
mol=loadmol2 {0}
solvateOct mol TIP3PBOX 14.0
set default PBradii {3}
saveAmberParm mol {1}/{2}-ligand.solv.top {1}/{2}-ligand.solv.crd
savepdb mol {1}/{2}-ligand.solv.amber.pdb'''.format( newligandfile, dir, ligand_name, radii))
        fhandle.close()

def add_belly(fhandle, protein_belly):
    fhandle.write('''\
  ibelly=1, bellymask = ":%s",
  /\n''' % protein_belly)
    fhandle.close()

def add_restraints(fhandle, ligand_restraints, protein_belly, restraint_k=0.5):
    if ligand_restraints==True and protein_belly!=None: 
        fhandle.write('''\
  ntr = 1, 
  restraintmask = ":MOL", restraint_wt = {0},
  ibelly=1, bellymask = ":{1}",
  /\n'''.format(restraint_k, protein_belly))
    elif ligand_restraints==True:
        fhandle.write('''\
  ntr = 1, 
  restraintmask = ":MOL", restraint_wt = {0},
  /\n'''.format(restraint_k))
    elif protein_belly!=None:
        fhandle.write('''\
  ibelly=1, bellymask = ":{0}",
  /\n'''.format(protein_belly))
    else:
        fhandle.write('''\
  ntr = 0,
  /\n''')     #to get around annoying issue with newline character  
    fhandle.close()

def write_simulation_input(md, dir, prefix, gb_model=None, gbmin=False, protein_belly=None, maxcycles=50000, steps=100000):
    filename='%s/%s.in' % (dir, prefix)
    if gbmin==False:
        fhandle=open(filename, 'w')
        fhandle.write('''\
minimization
  &cntrl
  imin = 1, maxcyc = %s, ntmin = 1,
  ncyc   = 100, 
  ntx = 1, ntc = 1, ntf = 1,
  ntb = 1, ntp = 0,
  ntwx = 1000, ntwe = 0, ntpr = 1000,
  drms   = 0.6,
  cut = 12.0,''' % maxcycles)
        if protein_belly!=None:
            add_belly(fhandle, protein_belly)
    if gbmin==True:
        fhandle=open(filename, 'w')
        fhandle.write('''\
minimization
  &cntrl
  imin = 1, maxcyc = %s, ntmin = 1,
  ncyc   = 100, 
  ntx = 1, ntc = 1, ntf = 1,
  ntb = 1, ntp = 0,
  ntwx = 1000, ntwe = 0, ntpr = 1000,
  igb=%s,
  drms   = 0.6,
  cut = 12.0,''' % (maxcycles, gb))
        if protein_belly!=None:
            add_belly(fhandle, protein_belly)
    if md==True:
        filename='%s/%s.in' % (dir, prefix)
        fhandle=open(filename, 'w')
        fhandle.write('''\
nvt equilibration with Langevin therm, SHAKE Hbonds
  &cntrl
  imin = 0, ntx = 1, irest = 0, nstlim = %s,
  ntt=2, temp0 = 298.15, tempi = 0, ig = -1,
  ntc = 2, ntf = 2, dt = 0.002,
  ntb = 1, ntp = 0, 
  ntwx = 1000, ntwe = 0, ntwr = 1000, ntpr = 1000,
  cut = 10.0, iwrap = 1,
  nscm = 100,''' % steps)
        if protein_belly!=None:
            add_belly(fhandle, protein_belly)

def write_mmgbsa_input(filename, model, start, interval, finish):
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
    if 'mol2' not in outconf:
        print "need mol2 complex-miminized ligand conf"
        sys.exit()
    fhandle=open(filename, 'w')
    fhandle.write('''\
trajin %s restart
strip !:MOL
trajout %s mol2''' % (inconf, outconf)) 
    fhandle.close()

