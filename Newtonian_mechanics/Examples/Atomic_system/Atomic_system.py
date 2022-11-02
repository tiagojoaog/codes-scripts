#!/usr/bin/env python3

from classical_mechanics import *

G = 1e-45  #gravitational constant in au
Body_list=[]

proton1  = kutta4(c=color.red,q=1,mass=1836.1526734,move=False,label='Proton 1',body_list=Body_list)
proton2  = kutta4(x=0.2,q=1,c=color.red,mass=1836.1526734,move=False,label="Proton 2",body_list=Body_list)
neutron1 = kutta4(x=0.1,y=0.1,q=0,c=color.white,mass=1836.1526734,move=False,label="Neutron 1",body_list=Body_list)
neutron2 = kutta4(x=0.1,y=-0.1,q=0,c=color.white,mass=1836.1526734,move=False,label="Neutron 2",body_list=Body_list)

electron1 = kutta4(x=1,r=10,q=-1,vz=1,c=color.yellow,label="Electron 1",body_list=Body_list) #r=1 au is 1 bohr = 0.58e-10m,v = 1au is 2.18e6 m/s = 0.7267% c (speed of light) 
electron2 = kutta4(x=2,r=10,q=-1,vz=0.05,c=color.yellow,label="Electron2",body_list=Body_list)
electron3 = kutta4(x=-2,r=10,q=-1,vz=-0.05,c=color.yellow,label="Electron3",body_list=Body_list)
electron4= kutta4(x=-1,r=10,q=-1,vz=-0.5,c=color.yellow,label="Electron 2",body_list=Body_list)
electron5= kutta4(y=2,r=10,q=-1,vy=0.03,c=color.yellow,label="Electron5",body_list=Body_list)

RM,M=center_of_mass(Body_list)

Atomic_system=orbit(G=G,Ke=1,frames=200,dt=0.001,RCM=RM,capture=False,is_adaptive=False,bodies=Body_list) #dt=1au= 1/41 fs (1/41 e-15 s)
Atomic_system.run(clear_cache=True)                                                                       #and t is 1/365.25 au = 1/(41*365.25) fs
Atomic_system.close()                                                                                     #this ratio is used to convert days to years
                                                                                                          #on celestial orbits, but is kept for atomic
                                                                                                          #systems as well. Ke=1 is in atomic units.
