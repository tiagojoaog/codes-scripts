#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
#from myconstraints import stretchcombo
from ase.constraints import FixAtoms, FixBondLength, FixedPlane
from ase.io import *
import sys
import os
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
##from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np


try:
  atoms=read('in.traj')
except:
  atoms=read('in.xyz')
  atoms.center(vacuum=9)


try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'


atoms.set_constraint()
atoms.pbc=True

calc = vasp_calculator.Vasp(encut=800,
                        xc='PBE',
                        gga='PE',
                        kpts  = (1,1,1),
                        ncore=4,
                        istart=1,
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        nelm=250,
                        algo = 'fast',
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-8,  #energy conv.
                        prec='Accurate',
                        nsw=1, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        lsol=True,
                        eb_k=80.0,
                        sigma_k=0.6,
                        nc_k=0.0025,       # converge these values to simulate the explicit solvation case
                        tau=0.000525,      # converge these values to simulate the explicit solvation case
                        ivdw=11,
                        ispin=1)
atoms.set_calculator(calc)
trajfile='final.traj'
dyn = BFGS(atoms, logfile=submitdir+'relax.log', trajectory=submitdir+'relax.traj')
dyn.run(fmax=0.01)
write(trajfile,atoms)
relaxed_atoms = read(trajfile)
e=relaxed_atoms.get_potential_energy()
print('final energy',e)
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
