#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
#from ase.structure import bulk
from ase.lattice.spacegroup import crystal
from ase.io import *
import sys
import os
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.optimize import BFGS
import numpy as np
from ARPESS.constraints_fixedcombo import stretchcombo
from ARPESS.arpess import arpess
#from ase_devel.optimize    import scan_constraint
from ase.constraints import FixAtoms, FixBondLength, FixedPlane

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

#els=atoms.get_chemical_symbols()
#npt=els.index('Cu')
#pole=[0. for ii in range(len(els))]
#pole[npt]=1.0
#atoms.set_initial_magnetic_moments(pole)

#atoms.set_constraint()

original_constraints=atoms.constraints

liste=[]
liste+=[[45,20,1.0]] # O1-C
liste+=[[21,20,-1.0]]
initial_value=0.0

# define constraint
cc=stretchcombo(initial_value, liste, atoms)
# adjust structure according to the intial value (in this example the intial value is not modified).
atoms.set_positions(cc.return_adjusted_positions()) 
# set all constrainsts, keeping old constraints (there are no old constraints in this example).  
atoms.set_constraint([cc])


calc = vasp_calculator.Vasp(encut=400.0,
                        xc='PBE',
                        gga='PE',
                        kpts  = (4,4,1),
                        ismear=0,
                        gamma=True,
                        nbands=320,
			ncore=16,
                        ivdw=11,
                        nelm=150,
                        algo = 'fast',
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-7,  #energy conv.
                        lreal=False,
                        prec='Accurate',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        ispin=2)     #TURN SPIN-POLARIZED CALCS ON

atoms.set_calculator(calc)

trajfile='final.traj'
dyn = BFGS(atoms, logfile=submitdir+'relax.log', trajectory=submitdir+'relax.traj')
dyn.run(fmax=0.01)
atoms.set_constraint(original_constraints)
write(trajfile,atoms)
relaxed_atoms = read(trajfile)
e=relaxed_atoms.get_potential_energy()
f=relaxed_atoms.get_forces(apply_constraint=True,md=False)
print("")
print('final energy',e)
print("")
print('forces:')
print(f)
fmax=[abs(max(i, key=abs)) for i in f]
fmax=float(max(fmax))
print("")
print('fmax',str(fmax))
f=open('final.e','w')
f.write('final energy (eV) '+str(e)+' fmax '+str(fmax)+'\n')
f.close()

#dyn = arpess(atoms, constraintlist=[cc],a=initial_value,maxstep=0.1,adaptive_threshold=None,linear_interpol=False,do_cubic=False, log_prefix=submitdir,traj_prefix=submitdir)
#dyn.run(fmax=0.01)

