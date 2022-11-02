#!/usr/bin/env python3

from classical_mechanics import *
import matplotlib.pyplot as plt

G = 6.67e-11  #gravitational constant in AU


Body_list1=[]

main_asteroid_rk4   = kutta4(body_list=Body_list1,r=35000,mass=5.4e11,draw=True,label="Asteroid 1",c=color.purple)
second_asteroid_rk4 = kutta4(x=1190,r=15000,vz=0.174,c=color.green,mass=5e9,label="Asteroid 2",body_list=Body_list1,draw=True)  # 1AU -distance from earth to sun     

RM,M=center_of_mass(Body_list1)


system1=orbit(G=G,Ke=0,frames=200,draw=True,dt=10,traj_show_save='save',RCM=RM,capture=False,is_adaptive=False,bodies=Body_list1) #dt=1 is 1 year
system1.run(timelimit=12*3600/365,clear_cache=True)  

Body_list1=[]

main_asteroid_rk4   = kutta4(body_list=Body_list1,r=35000,mass=5.4e11,draw=True,label="Asteroid 1",c=color.purple)
second_asteroid2_rk4 = kutta4(x=1190,r=15000,vz=0.174-2e-3,c=color.red,mass=5e9,label="Asteroid 2(after impact)",body_list=Body_list1,draw=True) 


system2=orbit(G=G,Ke=0,frames=200,draw=True,dt=10,traj_show_save='save',RCM=RM,capture=False,is_adaptive=False,bodies=Body_list1) #dt=1 is 1 year
system2.run(timelimit=12*3600/365,clear_cache=True)

fig=plt.figure()
ax=plt.axes(projection='3d')

for i in system1.properties:
    positions = system1.properties[i]['pos']
    X=[]
    Y=[]
    Z=[]
    for j in range(len(positions)):
        X.append(system1.properties[i]['pos'][j].x)
        Y.append(system1.properties[i]['pos'][j].y)
        Z.append(system1.properties[i]['pos'][j].z)
    ax.plot3D(X,Y,Z,label=i)


for i in system2.properties:
    positions = system2.properties[i]['pos']
    X=[]
    Y=[]
    Z=[]
    for j in range(len(positions)):
        X.append(system2.properties[i]['pos'][j].x)
        Y.append(system2.properties[i]['pos'][j].y)
        Z.append(system2.properties[i]['pos'][j].z)
    ax.plot3D(X,Y,Z,label=i)

ax.view_init(155,-95)
plt.legend(loc='best')
plt.show()


