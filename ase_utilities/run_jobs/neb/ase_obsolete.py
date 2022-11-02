#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from ase.io import read,Trajectory
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.optimize import BFGS
from matplotlib import pyplot as plt
import ase
import sys
import os
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.optimize import BFGS
import numpy as np
import subprocess

def swap_atoms(atoms0,swap):
  atoms=atoms0.copy()
  els=atoms.get_chemical_symbols()
  xyz=atoms.get_positions()
  news={}
  for ii in swap:
    news[ii]=[els[swap[ii]],list(xyz[swap[ii]])]
  for ii in swap:
    els[ii] = news[ii][0]
    xyz[ii] = news[ii][1]
  atoms.set_positions(xyz)
  atoms.set_chemical_symbols(els)
  return atoms

try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'

x=0

try:
  i=0
  images=[]
  trial_initial = read('0_custom.traj')
  #trial_initial.set_constraint()
  trig=1
  while trig==1:
    try:
      a = read(str(i)+'_custom.traj')
      #a.set_constraint()
      images.append(a)
      i+=1
    except:
      trig=0
  nimages=i
  x=1
  print("Custom path found with "+str(nimages)+" images.")
except:
  print("No custom path found.")

try:
  with open("out_restart.txt") as al:
    pass
  os.system("grep 'computing' out_restart.txt |awk '{print $2}'|tail -1>nimages")
  with open("nimages") as fil:
    for line in fil:
      nimages=int(line)
  images=ase.io.read('pre-neb.traj',index=slice(-nimages,None))
  ase.io.write('test.traj',images)
  x=1
except:
  print("pre-neb.traj or out_restart.txt doesnt exist.")

try:
  with open("out_restart.txt") as al:
    pass
  os.system("grep 'computing' out_restart.txt |awk '{print $2}'|tail -1>nimages")
  with open("nimages") as fil:
    for line in fil:
      nimages=int(line)
  images=ase.io.read('neb.traj',index=slice(-nimages,None))
  ase.io.write('test.traj',images)
  x=1
except:
  print("neb.traj or out_restart.txt doesnt exist.")  

if x==0:
  initial = read('initial.traj')
  final = read('final.traj')
  #initial.set_constraint()
  #final.set_constraint()

#swap={}
## final  initial position
#swap[113]=111
#swap[112]=110
#swap[111]=109
#swap[109]=112
#swap[110]=113

#final=swap_atoms(final,swap)
ispin=1
nupdown=-1  #not set (-1) 


if x==0:
  images= [initial]
  nimages=13

  for i in range(nimages-2):
        image = initial.copy()
        images.append(image)

  images.append(final)

subprocess.call("echo 'computing "+str(nimages)+" images'",shell=True,stdout=None,cwd=".")
images_complete = []
for i in range(len(images)):
  a=images[i]
  els=a.get_chemical_symbols()
  if ispin==2 and nupdown==1:
    try:
      npt=els.index('Cu')
      pole=[0. for ii in range(len(els))]
      pole[npt]=float(nupdown)
      a.set_initial_magnetic_moments(pole)
    except:
      print("Copper doesnt exist, magmom allocated somewhere else.")
  a.calc=vasp_calculator.Vasp(encut=530,
                        xc='PBE',
                        ivdw=12,
                        gga='PE',
                        kpts  = (4,4,1),
                        ncore=8,
                        nupdown=nupdown,
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        nelm=900,
                        algo = 'fast',
                        sigma = 0.1,
                        ibrion=-1,
                        ediffg=-0.01,  # forces
                        ediff=1e-7,  #energy conv.
                        prec='Accurate',
                        nsw=0, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=ispin)
  images_complete.append(a)


#neb = NEB(images_complete,climb=False,parallel=True)
neb = NEB(images_complete,climb=False,k=0.1)
#neb = NEB(images_complete,climb=True,method='improvedtangent')

idpp=False
idpp=True

if x==0:
  if idpp:
    neb.interpolate("idpp")
  else:
    neb.interpolate()
iii=0


runtrue=False
runtrue=True

traj = Trajectory('starting_trajectory.traj','w',neb.images[0])
for ii in neb.images:
  iii+=1
  #if iii==1:
  #  continue
  stri=str(iii)
  if len(stri)==1:
    stri='0'+stri
  traj.write(atoms=ii,mode='a')
traj.close()

if runtrue:
  qn = BFGS(neb, trajectory=submitdir+'pre-neb.traj',logfile=submitdir+'pre-neb.log',)
  qn.run(fmax=0.1)
  neb.climb=True
  qn = BFGS(neb, trajectory=submitdir+'neb.traj',logfile=submitdir+'neb.log',)
  qn.run(fmax=0.01)
  #neb.images[int(len(neb.images)/2)].write('ts.traj')
  os.system("echo 1>tmpone.txt")
  os.system("lookat_neb.py neb.traj<tmpone.txt")
  os.system("rm tmpone.txt")
