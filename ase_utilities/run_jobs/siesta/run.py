#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
from ase.io import *
import sys 
import os
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
from ase.calculators.dftd3 import DFTD3
import numpy as np
from ase.calculators.siesta.siesta import Siesta_Relax_with_ASE as Siesta
from ase.units import Ry

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



extra_arguments = {'DM.MixingWeight'               : 0.1,
                   'DM.Require.Energy.Convergence' : True,  
                   'DM.Energy.Tolerance'           : '10e-8 eV',    #Energy diff threshold in eV
                   'MaxSCFIterations'              : 300,           #SCF iterations
                   'MD.TypeOfRun'                  : 'CG',          #Defines md or relax simulations (CG,Broyden, FIRE).
                   'MD.MaxForceTol'                : '0.01 eV/Ang', #threshold for forces in eV/A (ignored if Forces is used)
                   'MD.VariableCell'               : False,         #to converge lattice parameters
                   'MD.NumCGsteps'                 : 0,             #ionic steps, 0 because relax,tion is done with ASE
                   'Diag.ParallelOverK'            :False,          #parallel over orbitals are usually preferred. If you have a small unit cell, with just a few atoms but much higher 
                                                                    #kpoint sampling, compared to number of processors (e.g. primitive Cu cell kpts 18x18x18) then set this key to True
                    }

calc = Siesta( label         = 'SIESTA_run',    #do not change this name, ase command depends on it
               xc            = 'PBE',           #exchange and correlation functional
               mesh_cutoff   = 200 * Ry,        #mesh cutoff in eV (200 Ry in this case)
               energy_shift  = 0.005 * Ry,       #energy shift in eV
               basis_set     = 'TZP2',          #basis set (in this case triple zeta with double polarization)
               kpts          = [4, 4, 1],       #kpoint sampling
               spin          = 'collinear', # “non-polarized”,”collinear”,"non-collinear" and "spin-orbit" are allowed
               fdf_arguments = extra_arguments)


calc = DFTD3(dft=calc,xc='PBE',damping='zero')          #Comment this for no D3 dispersion
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


