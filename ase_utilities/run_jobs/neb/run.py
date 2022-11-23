#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase,os,sys,subprocess,ase.io.vasp
from ase.io import *
from ase import Atoms
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
import numpy as np
from ase.neb import NEB

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


def sort_species(atom,symbol_count):
    Z       = atom.get_chemical_symbols()
    POS     = atom.get_positions()
    cell    = atom.get_cell()
    new_Z   = []
    new_POS = []

    if atom.constraints:
      for constr in atom.constraints:
        if isinstance(constr, FixAtoms):
          fixed = constr.todict()['kwargs']['indices']

    indices = []
    for sym, count in symbol_count:
        index1 = 0
        index2 = 0
        for z,pos in zip(Z,POS):
            if sym == z:
                if index1 in fixed:
                  indices.append(index2)
                new_Z.append(z)
                new_POS.append(pos)
                index2+=1
            index1+=1

    return Atoms(symbols=new_Z, positions=new_POS, cell=cell, constraint = FixAtoms(indices=indices))

def dict_to_list(_dict):
    _list=[]
    for name,value in _dict.items():
        _list.append([name,value])
    return _list

def species(atom):
    _dict   = {}
    for i in atom.get_chemical_symbols():
        if i not in _dict:
            _dict[i]=1
        else:
            _dict[i]+=1
    _list = dict_to_list(_dict)
    return _list

try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'


'''

SCRIPT STARTS HERE

'''

nimages=12  # total number of images (fixed initial and final ones also count)
idpp='idpp' # options are 'idpp' or ''
run=True    # if False, it will only generate the interpolated POSCAR files


POSCAR = True
for i in range(nimages):
    if os.path.exists(f'{i:0>2}/POSCAR'):
        pass
    else:
        POSCAR = False
      
if not POSCAR:
  initial = read('initial.traj')
  final   = read('final.traj')
  images= [initial]
  for i in range(nimages-2):
        image = initial.copy()
        images.append(image)

  images.append(final)

  neb = NEB(images,climb=True,k=0.1) #climb and k value are redundant here, ase NEB is only used for the interpolation of images, climb and k are defined inside the calculator
  if idpp:
    neb.interpolate(idpp)
  else:
    neb.interpolate()

  for index,image in enumerate(neb.images):
      os.system(f'mkdir {index:0>2}')
      _species = species(image)
      image = sort_species(image,_species)                                                          # sort_species is used to order positions according to _species list, which will be defined as symbol_count later
      ase.io.vasp.write_vasp(f'{index:0>2}/POSCAR',image,symbol_count=_species)                     # symbol count is used to write POSCARS in compact notation e.g. H C O rather than H C O C H (3 species only)
                                                                                                    # this is relevant because otherwise vasp will take the former example as 5 species instead of 3
                                                                                                    # and will crash because potcar only has 3 available species.

extra_string=''
if POSCAR:
  extra_string+='POSCARs read from ##/POSCAR files.'
else:
  extra_string+=f'POSCARs not read. Generating POSCARs from ASE NEB({idpp}) interpolation method.'

subprocess.call(f"echo 'computing {str(nimages)} images. {extra_string}'",shell=True,stdout=None,cwd=".")


if run:
   atoms=read('00/POSCAR')
   calc = vasp_calculator.Vasp(encut=400,
                        xc='PBE',
                        gga='PE',
                        ncore=8,
                        isif=2,
                        images=nimages-2,     #start NEB
                        spring=-5.0,
                        ichain=0,
                        lclimb=True,   #end NEB
                        ivdw=11,
                        kpts  = (1,1,1),
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        sigma = 0.1,
                        nelm=250,
                        algo = 'fast',
                        ibrion=1,    # -1 for no relaxation with vasp, 1 otherwise
                        ediffg=-0.01,  # forces
                        ediff=1e-6,  #energy conv.
                        prec='Accurate',
                        nsw=500, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=1)
   atoms.set_calculator(calc)
   e=atoms.get_potential_energy()
   print('final energy',e)
   f=open('final.e','w')
   f.write(str(e)+'\n')
   f.close()


