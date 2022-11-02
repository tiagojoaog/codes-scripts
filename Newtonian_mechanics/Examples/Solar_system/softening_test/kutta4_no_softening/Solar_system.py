#!/usr/bin/env python3

from classical_mechanics import *

G = 2.959122083e-4  # Gravitational constant in AU
softening = 0       # Scales down forces for very close objects - to avoid approaching infinity (1E-8 is ok for this system). 
                    # On this example, softening is only used on Sun and Earth, to prevent infinity for aproaching objects (to these two bodies).
Body_list=[] # Everytime an object is created (e.g. Sun), that object is appended to this list.

sun       = kutta4(softening=softening,body_list=Body_list,light=True,c=color.yellow,draw=False) # draw=False boosts the calculation, use this in every object to disable all animations
earth     = kutta4(x=1,r=10,softening=softening,texture=textures.earth,mass=3.003e-06,label="Earth",body_list=Body_list,draw=False)  # 1AU -distance from earth to sun, softening prevents rapid changes in acceleration (earth-moon system)     
jupiter   = kutta4(x=5.2,r=17,c=color.orange,mass=9.551e-04,label="Jupiter",body_list=Body_list,draw=False)
moon      = kutta4(x=1.00256954,r=10,c=color.white,mass=3.69369e-08,label="Moon",body_list=Body_list,draw=False)
venus     = kutta4(x=0.723,r=10,c=color.red,mass=3.01035735e-08,label="Venus",body_list=Body_list,draw=False)
asteroid  = kutta4(x=1.5,r=7,vz=0.005,c=color.green,mass=1e-08,label="asteroid",body_list=Body_list,draw=False) #velocities of these bodies are explicitly stated here.
asteroid2 = kutta4(y=2,r=10,vx=0.005,vz=0.005,c=color.purple,mass=8.4e-08,label="asteroid 2",body_list=Body_list,draw=False)

RM,M=center_of_mass(Body_list)

earth.circular_orbit(center=RM,bodies=Body_list,MASS=M,G=G)   # v=1 is 1 AU/year.
jupiter.circular_orbit(center=RM,bodies=Body_list,MASS=M,G=G) # For the solar system, circular orbit - velocity calculator accurately predicts the correct starting velocities of the bodies orbiting the sun.
moon.circular_orbit(body=earth,v=earth.v,G=G)# This calculates the moon's velocity with respect to earth, and adds up earth velocity to match up earth's velocity.
venus.circular_orbit(center=RM,bodies=Body_list,MASS=M,G=G)   # circular orbit - velocity calculator only works well on a 2D xz-plane.


Solar_system=orbit(G=G,Ke=0,frames=False,dt=0.05,RCM=RM,traj_show_save='save',capture=False,is_adaptive=True,bodies=Body_list,draw=False) #dt=1 is 1 day. If is_adaptive=True, dt is ignored and is
                                                                                                                                          #                calculated in each iteration  
Solar_system.run(threshold_adaptive=1E-15,max_dt=1,timelimit=2,clear_cache=True)  #Threshold_adaptive is the maximum error on conservation of total energy and angular momentum
Solar_system.close()                                                                #to calculate dt. max_dt is the maximum limit that the step can take. These variables belong to adaptive step.
                                                                                    #If timelimit is stated, the program will finish after the stated time (in years).
                                                                                    #If clear_cache True, the program will delete all the pkl files formed at the end of execution.
