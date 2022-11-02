#!/usr/bin/env python3

from classical_mechanics import *
import matplotlib.pyplot as plt

G = 2.959122083e-4  #gravitational constant in AU
Body_list1=[]

sun_euler   = euler(body_list=Body_list1,draw=False,c=color.yellow)
earth_euler = euler(x=1,r=10,c=color.blue,mass=3.003e-06,label="Earth",body_list=Body_list1,draw=False)  # 1AU -distance from earth to sun     
jupiter_euler = euler(x=5.2,r=17,c=color.orange,mass=9.551e-04,label="Jupiter",body_list=Body_list1,draw=False)
moon_euler  = euler(x=1.00256954,r=10,c=color.white,mass=3.69369e-08,label="Moon",body_list=Body_list1,draw=False)  
venus_euler = euler(x=0.723,r=10,c=color.red,mass=3.01035735e-08,label="Venus",body_list=Body_list1,draw=False)
asteroid_euler= euler(x=1.5,r=7,vz=0.005,c=color.green,mass=1e-08,label="asteroid",body_list=Body_list1,draw=False)
asteroid2_euler= euler(y=2,r=10,vx=0.005,vz=0.005,c=color.purple,mass=8.4e-08,label="asteroid 2",body_list=Body_list1,draw=False)

RM,M=center_of_mass(Body_list1)

earth_euler.circular_orbit(center=RM,bodies=Body_list1,MASS=M,G=G)      #v=1 is 1 AU/year
jupiter_euler.circular_orbit(center=RM,bodies=Body_list1,MASS=M,G=G)
moon_euler.circular_orbit(body=earth_euler,v=earth_euler.v,G=G)
venus_euler.circular_orbit(center=RM,bodies=Body_list1,MASS=M,G=G)

Solar_system1=orbit(G=G,Ke=0,frames=False,draw=False,dt=0.05,traj_show_save='save',RCM=RM,capture=False,is_adaptive=False,bodies=Body_list1) #dt=1 is 1 year
Solar_system1.run(timelimit=24.00,clear_cache=True)

Body_list2=[]

sun_rk2   = kutta2(body_list=Body_list2,draw=False,c=color.yellow)
earth_rk2 = kutta2(x=1,r=10,c=color.blue,mass=3.003e-06,label="Earth",body_list=Body_list2,draw=False)  # 1AU -distance from earth to sun     
jupiter_rk2 = kutta2(x=5.2,r=17,c=color.orange,mass=9.551e-04,label="Jupiter",body_list=Body_list2,draw=False)
moon_rk2  = kutta2(x=1.00256954,r=10,c=color.white,mass=3.69369e-08,label="Moon",body_list=Body_list2,draw=False)  
venus_rk2 = kutta2(x=0.723,r=10,c=color.red,mass=3.01035735e-08,label="Venus",body_list=Body_list2,draw=False)
asteroid_rk2= kutta2(x=1.5,r=7,vz=0.005,c=color.green,mass=1e-08,label="asteroid",body_list=Body_list2,draw=False)
asteroid2_rk2= kutta2(y=2,r=10,vx=0.005,vz=0.005,c=color.purple,mass=8.4e-08,label="asteroid 2",body_list=Body_list2,draw=False)

RM,M=center_of_mass(Body_list2)

earth_rk2.circular_orbit(center=RM,bodies=Body_list2,MASS=M,G=G)      #v=1 is 1 AU/year
jupiter_rk2.circular_orbit(center=RM,bodies=Body_list2,MASS=M,G=G)
moon_rk2.circular_orbit(body=earth_rk2,v=earth_rk2.v,G=G)
venus_rk2.circular_orbit(center=RM,bodies=Body_list2,MASS=M,G=G)

Solar_system2=orbit(G=G,Ke=0,frames=False,draw=False,dt=0.05,traj_show_save='save',RCM=RM,capture=False,is_adaptive=False,bodies=Body_list2) #dt=1 is 1 year
Solar_system2.run(timelimit=24.00,clear_cache=True)  

Body_list3=[]

sun_rk4   = kutta4(body_list=Body_list3,draw=False,c=color.yellow)
earth_rk4 = kutta4(x=1,r=10,c=color.blue,mass=3.003e-06,label="Earth",body_list=Body_list3,draw=False)  # 1AU -distance from earth to sun     
jupiter_rk4 = kutta4(x=5.2,r=17,c=color.orange,mass=9.551e-04,label="Jupiter",body_list=Body_list3,draw=False)
moon_rk4  = kutta4(x=1.00256954,r=10,c=color.white,mass=3.69369e-08,label="Moon",body_list=Body_list3,draw=False)
venus_rk4 = kutta4(x=0.723,r=10,c=color.red,mass=3.01035735e-08,label="Venus",body_list=Body_list3,draw=False)
asteroid_rk4= kutta4(x=1.5,r=7,vz=0.005,c=color.green,mass=1e-08,label="asteroid",body_list=Body_list3,draw=False)
asteroid2_rk4= kutta4(y=2,r=10,vx=0.005,vz=0.005,c=color.purple,mass=8.4e-08,label="asteroid 2",body_list=Body_list3,draw=False)

RM,M=center_of_mass(Body_list3)

earth_rk4.circular_orbit(center=RM,bodies=Body_list3,MASS=M,G=G)      #v=1 is 1 AU/year
jupiter_rk4.circular_orbit(center=RM,bodies=Body_list3,MASS=M,G=G)
moon_rk4.circular_orbit(body=earth_rk4,v=earth_rk4.v,G=G)
venus_rk4.circular_orbit(center=RM,bodies=Body_list3,MASS=M,G=G)

Solar_system3=orbit(G=G,Ke=0,frames=False,draw=False,dt=0.05,traj_show_save='save',RCM=RM,capture=False,is_adaptive=False,bodies=Body_list3) #dt=1 is 1 year
Solar_system3.run(timelimit=24.00,clear_cache=True)  

plt.plot(Solar_system1.step_list,Solar_system1.totE_list,label='E. Euler')
plt.plot(Solar_system2.step_list,Solar_system2.totE_list,label='E. RK2')
plt.plot(Solar_system3.step_list,Solar_system3.totE_list,label='E. RK4')
plt.legend(loc='best')
plt.savefig('mechanical_energy_integrators.pdf')
plt.close()

plt.plot(Solar_system1.step_list,Solar_system1.ang_mom_list,label='Ang. mom. Euler')
plt.plot(Solar_system2.step_list,Solar_system2.ang_mom_list,label='Ang. mom. RK2')
plt.plot(Solar_system3.step_list,Solar_system3.ang_mom_list,label='Ang. mom. RK4')
plt.legend(loc='best')
plt.savefig('angmom_integrators.pdf')
plt.close()

Solar_system3.close()
