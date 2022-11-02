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
from ase.constraints import FixAtoms, FixBondLength, FixedPlane
from ase.constraints import FixAtoms
from ase.io.trajectory import PickleTrajectory
from ase.vibrations import Vibrations
from ase.calculators.espresso import Espresso
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np

#atoms=read('in.cif')
#atoms=read('CONTCAR')

try:
  atoms=read('in.traj')
except:
  atoms=read('in.xyz')
  atoms.center(vacuum=9)
  atoms.pbc=True
try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'

atomic_symbols = atoms.get_chemical_symbols()


pseudopotentials = {}

for atomic_symbol in atomic_symbols:
  pseudopotentials[atomic_symbol] = atomic_symbol+'.UPF'

cell = atoms.get_cell()
pos = atoms.get_positions()

z_max=-100
for i in pos:
    if i[2] > z_max:
        z_max = i[2]
emaxpos = round(((z_max + cell[2][2])/2)/cell[2][2],2)

ecutwfc = 50

calc = Espresso(ecutwfc = ecutwfc,
                ecutrho = 10*ecutwfc,
                pseudopotentials = pseudopotentials,
                input_dft = 'PBE',
                vdw_corr = 'grimme-d3', #DFT-D3
                dftd3_version = 3,  #zero damping
                kpts  = (5,5,1),
                koffset = (0,0,0),
                occupations = 'smearing',
                smearing = 'fd',
                tprnfor = True, #calculate forces
                tstress = True, #calculate stress
                degauss = 0.007,
                electron_maxstep = 250,
                conv_thr = 1e-9,  #energy conv.
                calculation = 'scf',
                diagonalization = 'david',
                mixing_mode = 'plain',
                mixing_beta = 0.5,
                scf_must_converge = False,
                tefield = True, # saw-like potential simulating electric field
                dipfield = True, # dipole correction
                edir = 3, #direction of dipole correction
                emaxpos = emaxpos, #saw-like potential value
                nstep = 1, # don't use the espresso internal relaxation, only use ASE
                real_space = False,
                verbosity = 'high',
                nspin = 1)

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

