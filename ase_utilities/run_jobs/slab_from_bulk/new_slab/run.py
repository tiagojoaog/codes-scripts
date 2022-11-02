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
#if not (os.path.exists('vdw_kernel.bindat')):
#                os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np
from vasp_stuff import get_nelec

#atoms=read('in.cif')i
#atoms=read('CONTCAR')

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
                        xc='PBE',
                        gga='PE',
                        ncore=8,
                        ivdw=11,
                        nbands=int(get_nelec(atoms)*0.5*1.5),
                        amin=0.01,
                        kpts  = (2,2,1),
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        sigma = 0.1,
                        nelm=250,
                        algo = 'fast',
                        ibrion=-1,    # -1 for no relaxation with vasp, 1 otherwise
                        ediffg=-0.005,  # forces
                        ediff=1e-7,  #energy conv.
#                        isif=3,
                        prec='Accurate',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=1)
                   #     idipol=3)
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
