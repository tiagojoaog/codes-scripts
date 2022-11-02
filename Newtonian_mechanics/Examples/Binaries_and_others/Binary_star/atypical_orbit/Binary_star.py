#!/usr/bin/env python3

from classical_mechanics import *


G = 2.959122083e-4  #gravitational constant in AU
Body_list=[]


Star1   = kutta4(body_list=Body_list,texture=textures.earth,label='Star 1')  # This creates a earth texture on this body. Any texture can be used as long as the file is available (not tested).
Star2   = kutta4(x=1,r=15,c=color.yellow,mass=0.7,light=True,label="Star 2",body_list=Body_list)  # light=True enables light effects from this object. This will be the light source.

RM,M=center_of_mass(Body_list)

Star1.circular_orbit(center=RM,bodies=Body_list,MASS=M,G=G) 
Star2.circular_orbit(center=RM,bodies=Body_list,MASS=M,G=G)  

Star1.v /= 2
Star2.v /= 2   #Cuts both calculated initial velocities (using circular_orbit definition) in half

binary=orbit(G=G,Ke=0,frames=200,dt=0.1,capture=False,is_adaptive=True,bodies=Body_list) #dt=1 is 1 year
binary.run(threshold_adaptive=1E-11,max_dt=1,clear_cache=True)
binary.close()
