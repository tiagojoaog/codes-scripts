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
import os,sys
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
                        gga='PE',
                        ivdw=11,
                        ldau=True,  #activate LDA+U
                        ldautype=2,
                        ldau_luj={'Fe':{'L':2,  'U':4.3, 'J':0.0},   #U=4.3 works well for alpha-Fe2O3 surfaces
                                 'O':{'L':-1, 'U':0.0, 'J':0.0}},
                        ldauprint=2,
                        lmaxmix=4,    # LDA+U  4 for up to d orbitals, 6 for up to f orbitals
                        lorbit=11,   
                        gamma=True,
                        kpts  = (4,4,1),
                        ismear=0,
                        nelm=950,
                        nbands = 320,
                        algo = 'fast',
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-5,  #energy conv.
                        lreal=False,
                        prec='Accurate',
                        nsw=1, # don't use the VASP internal relaxation, only use ASE
                        ispin=2)

atoms.set_calculator(calc)
e=atoms.get_potential_energy()
print('final energy',e)
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
