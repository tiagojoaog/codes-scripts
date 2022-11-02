#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
from ase.constraints import FixAtoms, FixBondLength, FixedPlane
from ase.io import *
import sys
import os
from ase.constraints import FixAtoms
from ase.io.trajectory import PickleTrajectory
from ase.vibrations import Vibrations
from ase.calculators.orca import ORCA
##from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np
from vasp_stuff import get_nelec


try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'

atoms=read('in.traj')
#atoms=read('CONTCAR')
#cell=atoms.get_cell()
#cell[0][0]=12.7
#cell[1][1]=12.7
#cell[2][2]=26.6
#atoms.set_cell(cell,scale_atoms=True)

params  = 'SCFConvForced SCFconv9 B3LYP Grid5 FinalGrid6 '  #replace grid and BLYP with DLPNO-CCSD(T) for CCSD(T) calcs
params += 'def2-TZVPP def2/JK RIJK D3ZERO KDIIS SOSCF Angs'     #add frozencore for HF/postHF, auxBasis is /C and /JK for hybrid funcs

'''pal'''
blocks  = '%pal   nprocs 32   end \n' 
''' '''

'''mdci'''
#blocks += '%MDCI  MaxCore 1000   end \n',    #for CCSD(T), replace MDCI with MP2 for MP2 calcs
''' '''

'''METHOD'''
blocks += '%METHOD   METHOD DFT end \n' #we need gradient, dont use runtyp SP
''' '''

'''scf'''
blocks += '%scf   HFTYP RKS \n  MaxCore 1000 \n DIISMaxEq 20 \n directresetfreq 1 \n' #replace RKS with RHF for H
blocks += 'SOSCFStart 0.0005 \n MAXITER 300 \n GUESS HUECKEL  end '  # SOSCFMaxIt x forces SOSCF to start after x iterations, disregarding SOSCFStart
''' '''

calc = ORCA(label='orcacalc',
       orcasimpleinput  = params,     
       orcablocks       = blocks,
       charge = 5, mult = 1,  
       )


atoms.set_calculator(calc)
e=atoms.get_potential_energy()
print('final energy',e)
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
