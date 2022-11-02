#!/usr/bin/env python3

import ase.io
import os,sys
import pickle
import numpy as np
import matplotlib.pyplot as plt


def all_combinations(numbers):
    return np.array(np.meshgrid(numbers, numbers, numbers)).T.reshape(-1,3)
    #combinations=[]
    #for i in numbers:
    #    for j in numbers:
    #        for k in numbers:
    #            if [i,j,k] not in combinations:
    #                combinations.append([i,j,k])
    #return np.array(combinations)

def shift_atom(cell,i,j,threshold):
    configurations=all_combinations([1,0,-1])
    modifying_matrix=np.dot(configurations,cell)

    for a in modifying_matrix:
        new_j=j+a
        if np.linalg.norm(new_j-i) < threshold:
            return new_j

def print_out(message,txt_file):
  print(message)
  neb.write(f'{message}\n')

def read_energy_from_file(filename):
  with open(filename) as e:
      return float([i for i in e][0])

def read_energy(efile,energy):
  if os.path.exists(efile):
      return read_energy_from_file(efile)
  else:
      return energy

def write_energy_to_file(filename,energy):
    e=open(filename,'w')
    e.write(energy)
    e.close()

def write_energy(efile,energy):
    if not os.path.exists(efile):
      write_energy_to_file(efile,str(energy))


just_plot = False



''' SCRIPT STARTS HERE '''

initial_energy  = False
final_energy    = False

for argument in sys.argv:
    if 'initial' in argument.split('='):
        initial_energy = float(argument.split('=')[1])
    elif 'final' in argument.split('='):
        final_energy = float(argument.split('=')[1])
    elif 'norun'==argument:
        just_plot=True
    elif 'h'==argument.lower() or 'help'==argument.lower():
        print('''Arguments:
                *initial=number   defines energy of the first image;
                *final=number     defines energy of the last image;
                *norun            prevents the script from running, if there is saved data e.g. neb.traj or neb.pkl, it will read the trajectory and plot the graphs (commands above will be ignored);
                *help             to check available commands.
                
                Note: by assigning initial and final once, energies will be stored, which will be read on later executions of this script.''')
        sys.exit(1)



if not just_plot:
    assume_energies = False
    neb             = open('neb.txt','w')
    

    if not initial_energy or not final_energy:
        assume_energies=True
        print_out('''
    
WARNING: Energies from initial and final images not defined. Script will assume adjacent image energies to be the same or will read from 'final.e' in 
         their folders if they exist. To define energies, either edit initial_energy=False and final_energy=False to their respective numbers in the script, 
         or execute the script as 'neb_vasp_read.py initial=number final=number'. By doing this, final.e will be stored in the first and last folders, 
         so that they can be read later.

         ''',neb)
    
    
    
    nimages        = 0
    run            = True
    
    iterations          =[]
    properties          ={}
    properties['E']     ={}
    properties['dE']    ={}
    properties['image'] ={}
    
    print_out('If iterations are higher than 200 per image, it may take more than 5 minutes. Be patient... \n\n\n',neb)

    while run:
        folder=f"{nimages:0>2}"
        if os.path.exists(folder):
            if not os.path.exists(f"{folder}/XDATCAR"):
                pass
            else:
                properties['E'][nimages]     ={}
                properties['dE'][nimages]    ={}
                properties['image'][nimages] ={}
    
                f=open(f"{folder}/stdout")
    
                for line in f:
                    _line = line.split(' ')
                    _line = list(filter(None,_line))
    
                    if 'F=' in _line:
                        niterations  =int(_line[0])
                        E            =float(_line[2])
                        dE           =float(_line[-1].split('=')[1])
    
                        template = ase.io.read(f'{folder}/POSCAR')
                        atoms    = ase.io.read(f'{folder}/XDATCAR',index=niterations-1)
                        cell     = template.get_cell()
    
                        threshold= min([np.linalg.norm(cell[0]),np.linalg.norm(cell[1]),np.linalg.norm(cell[2])])/2
                        
                        new_positions=[]
                        for i,j in zip(template.get_positions(),atoms.get_positions()):
                          distance=np.linalg.norm(j-i)
                          if abs(distance)>threshold:
                              j=shift_atom(cell,i,j,threshold)
                          new_positions.append(j)
                        
                        atoms.set_positions(new_positions)
    
                        properties['E'][nimages][niterations]            = E
                        properties['dE'][nimages][niterations]           = dE
                        properties['image'][nimages][niterations]        = atoms
    
                        iterations.append(niterations)
            nimages+=1
        else:
            run=False
    
    
    
    print_out(f'Number of images : {nimages}\n\n',neb)
    traj = ase.io.Trajectory('neb.traj','w')
    plot_energy = []
    plot_iter   = []
    
    initiale_file = '00/final.e'
    finale_file   = f'{nimages-1:0>2}/final.e'
    
    if not assume_energies:
        write_energy(initiale_file, initial_energy)
        write_energy(finale_file, final_energy)
    
    for n_iter in range(max(iterations)):
        n_iter+=1
        
        initial_energy=read_energy(initiale_file,properties['E'][1][n_iter])
        final_energy  =read_energy(finale_file,properties['E'][nimages-2][n_iter])
    
        properties['E'][0]                = {n_iter:initial_energy}
        properties['dE'][0]               = {n_iter:0}
        properties['image'][0]            = {n_iter:ase.io.read(f"00/POSCAR")}
        properties['E'][nimages-1]        = {n_iter:final_energy}
        properties['dE'][nimages-1]       = {n_iter:0}
        properties['image'][nimages-1]    = {n_iter:ase.io.read(f"{nimages-1:0>2}/POSCAR")}
        print_out(f'Iteration {n_iter:>4}:',neb)
        for n_im in range(nimages):
            delta=properties['dE'][n_im][n_iter]
            if abs(delta)>1E-4 or delta==0:
                delta_str=f'{delta:>10.4f}'
            else:
                delta_str=f'{delta:>10.2e}'
    
            print_out(f"   * image {n_im:0>2}: E = {properties['E'][n_im][n_iter]:>8.2f} eV   dE = {delta_str} eV",neb)
            plot_energy.append(properties['E'][n_im][n_iter])
            plot_iter.append(1+n_im+nimages*(n_iter-1))
            traj.write(atoms=properties['image'][n_im][n_iter],mode='a')
            if n_iter==max(iterations):
                properties['image'][n_im][n_iter].write(f'{n_im:0>2}/POSCAR')
        print_out('\n',neb)
   
    with open('neb.pkl','wb') as pkl:
      pickle.dump(plot_iter,pkl)
      pickle.dump(plot_energy,pkl)
    print_out('''
    
    POSCARS in folders have been updated. You can now restart the calculation.
    Creating energy plots (neb.pdf) and trajectories (neb.traj) ...
    
    ''',neb)

try:
  with open('neb.pkl','rb') as pkl:  
    plot_iter= pickle.load(pkl)
    plot_energy= pickle.load(pkl)
except:
  print("ERROR: Could not open pickle file, it can be damaged or does not exist at all.")
  sys.exit(10)

if just_plot:
    os.system("cat neb.txt")
plt.plot(plot_iter,plot_energy,marker='.')
plt.xlabel('Total Iterations')
plt.ylabel('Energy (eV)')
plt.savefig('neb.pdf')
os.system('ase gui neb.traj &')
plt.show()
plt.close()

