#!/usr/bin/env python3

from classical_mechanics import *
import matplotlib.pyplot as plt

G = 2.959122083e-4  #gravitational constant in AU


Body_list1=[]

sun_rk4   = kutta4(body_list=Body_list1,draw=False,c=color.yellow)
earth_rk4 = kutta4(x=1,r=10,c=color.blue,mass=3.003e-06,label="Earth",body_list=Body_list1,draw=False)  # 1AU -distance from earth to sun     
jupiter_rk4 = kutta4(x=5.2,r=17,c=color.orange,mass=9.551e-04,label="Jupiter",body_list=Body_list1,draw=False)
moon_rk4  = kutta4(x=1.00256954,r=10,c=color.white,mass=3.69369e-08,label="Moon",body_list=Body_list1,draw=False)  
venus_rk4 = kutta4(x=0.723,r=10,c=color.red,mass=3.01035735e-08,label="Venus",body_list=Body_list1,draw=False)
asteroid_rk4= kutta4(x=1.5,r=7,vz=0.005,c=color.green,mass=1e-08,label="asteroid",body_list=Body_list1,draw=False)
asteroid2_rk4= kutta4(y=2,r=10,vx=0.005,vz=0.005,c=color.purple,mass=8.4e-08,label="asteroid 2",body_list=Body_list1,draw=False)

RM,M=center_of_mass(Body_list1)

earth_rk4.circular_orbit(center=RM,bodies=Body_list1,MASS=M,G=G)      #v=1 is 1 AU/year
jupiter_rk4.circular_orbit(center=RM,bodies=Body_list1,MASS=M,G=G)
moon_rk4.circular_orbit(body=earth_rk4,v=earth_rk4.v,G=G)
venus_rk4.circular_orbit(center=RM,bodies=Body_list1,MASS=M,G=G)

Solar_system1=orbit(G=G,Ke=0,frames=False,draw=False,dt=0.05,traj_show_save='save',RCM=RM,capture=False,is_adaptive=False,bodies=Body_list1) #dt=1 is 1 year
Solar_system1.run(timelimit=24.00,clear_cache=True)  

Body_list2=[]

sun_rk4   = kutta4(body_list=Body_list2,draw=False,c=color.yellow)
earth_rk4 = kutta4(x=1,r=10,c=color.blue,mass=3.003e-06,label="Earth",body_list=Body_list2,draw=False)  # 1AU -distance from earth to sun     
jupiter_rk4 = kutta4(x=5.2,r=17,c=color.orange,mass=9.551e-04,label="Jupiter",body_list=Body_list2,draw=False)
moon_rk4  = kutta4(x=1.00256954,r=10,c=color.white,mass=3.69369e-08,label="Moon",body_list=Body_list2,draw=False)
venus_rk4 = kutta4(x=0.723,r=10,c=color.red,mass=3.01035735e-08,label="Venus",body_list=Body_list2,draw=False)
asteroid_rk4= kutta4(x=1.5,r=7,vz=0.005,c=color.green,mass=1e-08,label="asteroid",body_list=Body_list2,draw=False)
asteroid2_rk4= kutta4(y=2,r=10,vx=0.005,vz=0.005,c=color.purple,mass=8.4e-08,label="asteroid 2",body_list=Body_list2,draw=False)

RM,M=center_of_mass(Body_list2)

earth_rk4.circular_orbit(center=RM,bodies=Body_list2,MASS=M,G=G)      #v=1 is 1 AU/year
jupiter_rk4.circular_orbit(center=RM,bodies=Body_list2,MASS=M,G=G)
moon_rk4.circular_orbit(body=earth_rk4,v=earth_rk4.v,G=G)
venus_rk4.circular_orbit(center=RM,bodies=Body_list2,MASS=M,G=G)

Solar_system2=orbit(G=G,Ke=0,frames=False,draw=False,dt=0.05,traj_show_save='save',RCM=RM,capture=False,is_adaptive=True,bodies=Body_list2) #dt=1 is 1 year
Solar_system2.run(threshold_adaptive=1E-15,max_dt=0.15,timelimit=24.00,clear_cache=True)  

plt.plot(Solar_system1.step_list,Solar_system1.totE_list,label='Fixed dt (0.05)')
plt.plot(Solar_system2.step_list,Solar_system2.totE_list,label='Adaptive dt (1E-15 angmom and totE error)')
plt.legend(loc='best')
plt.savefig('mechanical_energy_adaptive.pdf')
plt.close()

plt.plot(Solar_system1.step_list,Solar_system1.ang_mom_list,label='Fixed dt (0.05)')
plt.plot(Solar_system2.step_list,Solar_system2.ang_mom_list,label='Adaptive dt (1E-15 angmom and totE error)')
plt.legend(loc='best')
plt.savefig('angmom_adaptive.pdf')
plt.close()

Solar_system2.close()
