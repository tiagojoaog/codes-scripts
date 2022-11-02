#! /usr/bin/env python
import ase
from ase.constraints import FixAtoms
#from ase.structure import bulk
from ase.lattice.spacegroup import crystal
from ase.io import *
import sys
import os
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
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


calc = vasp_calculator.Vasp(encut=400.0,
                        xc='PBE',
                        gga='PE',
                        kpts  = (1,1,1),
                        ismear=0,
			ncore=8,
                        ivdw=11,
                        nelm=250,
                        algo = 'fast',
                        laechg = True,
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-8,  #energy conv.
                        lreal='Auto',
                        prec='accurate',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        ispin=1)

atoms.set_calculator(calc)
e=atoms.get_potential_energy()
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
