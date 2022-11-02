#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
from ase.constraints import FixAtoms
#from ase.structure import bulk
from ase.lattice.spacegroup import crystal
from ase.io import *
from sys import path
import os
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np

if os.path.exists('in.traj'):
  atoms=read('in.traj')
else:
  atoms=read('in.xyz')
  atoms.center(vacuum=9)

calc = vasp_calculator.Vasp(encut=400.0,
                        xc='PBE',
                        gga='BF',
                        luse_vdw=True,
                        zab_vdw=-1.8867,
                        gamma=True,
                        kpts  = (1,1,1),
                        ismear=0,
                        nelm=250,
                        algo = 'fast',
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-8,  #energy conv.
                        lreal='Auto',
                        prec='normal',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        ispin=1)

atoms.set_calculator(calc)
e=atoms.get_potential_energy()
print('final energy',e)
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
