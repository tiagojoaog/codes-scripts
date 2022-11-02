#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
from myconstraints import stretchcombo
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
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'


atoms=read('in.traj')
atoms.center(vacuum=8,axis=2)

CM=atoms.get_center_of_mass(scaled=True) #scaled True to have lattice coordinates (e.g. 0.5a 0.5a 0.5a)
atoms.set_constraint()

calc = vasp_calculator.Vasp(encut=400,
                        xc='PBE',
                        ivdw=11,
                        gga='PE',
                        kpts  = (6,6,1),
			dipol = (CM[0],CM[1],CM[2]),
			ldipol = True,
                        idipol = 3,
                        ncore=8,
                                        # leave commented out for HSE calculations ... must be default
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        #nupdown=0.0,
                        nelm=250,
                        algo = 'fast',
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-8,  #energy conv.
                        prec='normal',
                        nsw=50, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=1)
atoms.set_calculator(calc)
trajfile='final.traj'
dyn = BFGS(atoms, logfile=submitdir+'relax.log', trajectory=submitdir+'relax.traj')
dyn.run(fmax=0.01)
write(trajfile,atoms)
relaxed_atoms = read(trajfile)
e=relaxed_atoms.get_potential_energy()
print 'final energy',e
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
