#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
from ase.constraints import FixAtoms, FixBondLength, FixedPlane
from ase.io import *
import sys
import os
from ase.constraints import FixAtoms
from ase.io.trajectory import PickleTrajectory
from ase.vibrations import Vibrations
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
##from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np
from vasp_stuff import get_nelec

print('start')
try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'

atoms=read('in.traj')
#atoms=read('CONTCAR')
#cell=atoms.get_cell()
#cell[0][0]=12.7
#cell[1][1]=12.7
#cell[2][2]=26.6
#atoms.set_cell(cell,scale_atoms=True)

calc = vasp_calculator.Vasp(encut=400,
                        xc='PBE',
                        #setups={'O': '_s'},
                        gga='PE',
                        kpts  = (4,4,1),
                        #kpar=1, # use this if you run on one node (most calculations).  see suncat confluence page for optimal setting
                        ncore=4,
                                        # leave commented out for HSE calculations ... must be default
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        #nupdown=0.0,
                        nelm=250,
			ivdw=11,
                        algo = 'fast',
                        #ismear = 1,
                        sigma = 0.1,
                        #algo = 'fast',
                        #ismear=-5,
                        #sigma=0.1,
                        ibrion=-1,
                        #ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-10,  #energy conv.
                        #nbands=100,
                        #nedos=2001,
                       # prec='Normal',
                        prec='Accurate',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        #lvtot=False,
                        lreal='Auto',
                 ##       isif=3, ## CELL EVERYTHING
                        ispin=1)
atoms.set_calculator(calc)
e=atoms.get_potential_energy()
print('final energy',e)
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
f=atoms.get_forces(apply_constraint=True)
print('norm of force',np.linalg.norm(f),', fmax_i = ',max(np.linalg.norm(fi) for fi in f))
print('>>> BEGIN print full force:')
for fi in f:
  print(fi[0],fi[1],fi[2])
print('<<< END print full force:')
if atoms.constraints:
  for constr in atoms.constraints:
    if isinstance(constr, FixAtoms):
      fixed = constr.todict()['kwargs']['indices']
else:
  fixed=[]
vibindices = [ ii for ii in range(len(atoms.get_positions())) if ii not in fixed]
vib = Vibrations(atoms,nfree=4,name=submitdir+'vib',delta=0.015,indices=vibindices)
vib.run()
vibenergies=vib.get_energies()
vib.summary(log='vib.txt')
for mode in range(len(vibindices)*3):
    vib.write_mode(mode)
