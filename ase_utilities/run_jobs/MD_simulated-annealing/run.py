#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
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

calc = vasp_calculator.Vasp(encut=400,
                        setups={'O': '_s'},
                        xc='PBE',
                        gga='PE',
                        ncore=8,
                        kpts  = (1,1,1),
                        ismear=0,
                        nelm=250,
                        ivdw = 11,
                        sigma = 0.1,
                        ediffg=-0.01,  # forces
                        isym=0,
                        nelmin=4,
                        algo = 'Fast',
######## MD settings
                        ibrion=0,
                        maxmix=40,
                        potim=1.0, # fs
                        nsw=int(1572),  # ionic steps for md T=TEBEG+(TEEND-TEBEG)×NSTEP/NSW
                        tebeg=524.0, # starting T
                        teend=0.0,    # final temperatures for simulated annealing, comment for standard MD
                        smass=-1.0, ## required for simulated annealing, 0 for standard MD with the Nose algorithm
                     #   nblock=50, #?
######## MD settings
                        lwave=False,
                        lcharg=False,
                        ediff=1e-7,  #energy conv.
                        lreal='Auto',
                        prec='normal',
                        ispin=1)

atoms.set_calculator(calc)
e=atoms.get_potential_energy()
