#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
#from ase.constraints import FixAtoms
#from ase.structure import bulk
#from ase.lattice.spacegroup import crystal
from ase.io import *
import sys
import os
import numpy as np


atoms = read('in.traj')

initial  = 3.5 
final    = 7.5
interval = 0.5

a=np.linspace(3.5,7.5,int((final-initial)/interval)+1)

here = os.getcwd()

for i in a:
  os.system(f'mkdir a_{i}')
  os.chdir(f'a_{i}')
  os.system('cp ../{in.traj,run.py,vasp.sh} .')
  atoms.set_cell(((i,0,0),(np.cos(np.pi*120/180)*i,np.sin(np.pi*120/180)*i,0),(0,0,28)))
  atoms.center(vacuum=20,axis=2)
  atoms.write('in.traj')
  os.system('sbatch vasp.sh')
  os.chdir(here)

