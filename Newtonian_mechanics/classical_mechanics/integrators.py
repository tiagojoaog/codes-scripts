#!/usr/bin/env python3

from vpython import *
from matplotlib import pyplot as plt
import numpy as np
import sys,os

__all__=["dummy","color","textures","euler","kutta2","kutta4"]

class dummy():
    def __init__(self):
        pos=vector(0,0,0)
class euler():
    def __init__(self,x=0,y=0,z=0,vx=0,vy=0,vz=0,softening=0,texture=None,append=True,light=None,v=None,pos=None,r=20,q=0,c=color.white,mass=1,label='Sun',draw=True,move=True,body_list=[]):
        self.properties          = {}    # __init__ initializes the object with all the properties. If properties not stated, default values are used. As shown above.
        self.images              = []
        self.mass                = mass
        self.label               = label
        self.light               = light
        self.move                = move
        self.q                   = q
        self.x                   = x
        self.y                   = y
        self.z                   = z
        self.vx                  = vx
        self.vy                  = vy
        self.vz                  = vz
        self.v                   = v
        self.pos                 = pos
        self.order               = 1
        self.e                   = softening

        if draw:
            self.object   = sphere(pos=self.pos, radius=r*4.65047e-3, color=c, texture=texture,
                                   make_trail=True, trail_type='points', interval=10, retain=50)
            if self.light:
                self.lamp = local_light(pos=self.pos,
                                                color=c)
        else:
            self.object     = dummy()    # Dummy is used in case we choose not perform animation
            self.lamp       = dummy()
            self.object.pos = self.pos
            self.lamp.pos   = self.pos

        self.properties['label']   = self.label  #object's properties
        self.properties['pos']     = self.pos
        self.properties['v']       = self.v
        self.properties['q']       = self.q
        self.properties['mass']    = self.mass
        self.properties['draw']    = draw
        self.properties['color']   = c
        self.properties['radius']  = r
        self.properties['light']   = self.light
        self.properties['texture'] = texture
        self.properties['move']    = self.move
        
        self.save_properties = self.properties.copy()
        self.images.append(self.save_properties)
        if append:
            body_list.append(self) #appends the entire object to a list. This list is then used to iterate over all objects to calculate accelerations and such.

    @property
    def v(self):
        return self._v
    @v.setter
    def v(self,v):
        if not v:
            self._v = vector(self.vx,self.vy,self.vz)
        else:
            self._v = v
    @v.deleter
    def v(self):
        self._v = None
    @property
    def pos(self):
        return self._pos
    @pos.setter
    def pos(self,pos):
        if not pos:
            self._pos = vector(self.x,self.y,self.z)
        else:
            self._pos = vector(pos)
    @pos.deleter
    def pos(self):                   
        self._pos = None
    def acceleration(self,b_object=None,G=0,Ke=0):    #acceleration definition
        r=b_object.pos-self.object.pos
        self.a  += (G * b_object.mass * r.hat / np.sqrt((mag(r)**4) +(self.e**2)))  - (Ke * self.q * b_object.q * r.hat / np.sqrt(((mag(r)**4) +(self.e**2))*self.mass))

    def update_vel_euler(self,dt=0):                   #updates velocities using euler method, as long as object is not constrained
        if self.move:
            self.v0 = self.v
            self.v += self.a*dt
    def update_pos(self,dt=0):                         #saves object's positions, updates animated object position as well as light source, if true
        if self.move:
            self.object.pos += ((self.v + self.v0)/2) * dt
            self.pos        = self.object.pos
            if self.light:
                self.lamp.pos = self.pos
        self.properties['pos'] = self.pos
        self.properties['v'] = self.v
        self.store()
    def store(self):                                  #saves properties to a list
        self.save_properties = self.properties.copy()
        self.images.append(self.save_properties)  
    def clear(self): # deletes the object's properties. This is called when these properties are flushed to pickle files (aka cache).
        del self.save_properties
        self.images = []
    def copy(self,label=None): #copies the entire object and its properties, used to calculate adaptive steps and kutta coefficients without updating properties of original bodies
        body = self.__class__(x=self.object.pos.x,y=self.object.pos.y,z=self.object.pos.z,v=self.v,q=self.q,mass=self.mass,label=label,move=self.move,draw=False,append=False)
        return body
    def energies(self,bodies,G=0,Ke=0,RCM=vector(1,0,0)): #calculates coulomb, gravitational potential energy, kinetic energy, angular momentum
        self.KE = (1./2)*self.mass*(mag(self.v))**2
        self.ANG_MOM=vector(0,0,0)
        self.GPE=0
        self.CPE=0
        self.force=vector(0,0,0)
        for i in bodies:
            if i!=self:
                r=i.object.pos-self.object.pos
                self.GPE   += -G*self.mass*i.mass/mag(r)
                self.CPE   += Ke*self.q*i.q/mag(r)
                self.force += (G * self.mass *i.mass * r.hat / np.sqrt((mag(r)**4) +(self.e**2)))  - (Ke * self.q * i.q * r.hat / np.sqrt(((mag(r)**4) +(self.e**2))))
                r_angular_mom=RCM-self.object.pos
                self.ANG_MOM+=cross(r_angular_mom,self.v*self.mass)
        self.PE=self.GPE+self.CPE
    def circular_orbit(self,center=None,MASS=None,body=None,bodies=None,v=None,G=0,Ke=0): #predicts initial velocities for circular orbits
        if center:
          m=MASS-self.mass        
          r=self.object.pos-center
          self.energies(bodies,G=G,Ke=Ke,RCM=center)
          circular_v = np.sqrt(mag(r)*mag(self.force)/self.mass)                               
          self.v += cross(vector(0,-circular_v,0),r.hat)  
        elif body:
          m=body.mass
          r=self.object.pos-body.object.pos
          circular_v = np.sqrt(1*m*G/mag(r))             #velocity of a circular orbit and escape velocity differ by ve=vc*sqrt(2)
          self.v += cross(vector(0,-circular_v,0),r.hat) 
        if v:
          self.v += v
    def __repr__(self):  # Outputs this when printing the object
        _list_=[]
        for i in self.images[-1]:
            _list_.append(i+':'+str(self.images[-1][i]))
        return '{0}({1})'.format(self.__class__.__name__, ', '.join(_list_))

class kutta2(euler): #inherits euler and adds or replaces definitions
    def __init__(self,x=0,y=0,z=0,vx=0,vy=0,vz=0,softening=0,append=True,light=None,v=None,pos=None,r=20,q=0,c=color.white,texture=None,mass=1,label='Sun',draw=True,move=True,body_list=[]):
        super().__init__(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,v=v,append=append,light=light,pos=pos,softening=softening,r=r,q=q,c=c,texture=texture,mass=mass,label=label,draw=draw,move=move,body_list=body_list)
        self.order=2
    def update_vel(self,params_vel_pos=None,dt=0): #kutta2 velocities
        if self.move:
            self.v0 = self.v
            k1=params_vel_pos[0]['copy_'+self.label]['a']
            k2=params_vel_pos[1]['copy_'+self.label]['a']
            self.v  += (k1+k2)*dt/2
class kutta4(euler):
    def __init__(self,x=0,y=0,z=0,vx=0,vy=0,vz=0,softening=0,append=True,light=None,v=None,pos=None,r=20,q=0,c=color.white,texture=None,mass=1,label='Sun',draw=True,move=True,body_list=[]):
        super().__init__(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,v=v,append=append,light=light,pos=pos,softening=softening,r=r,q=q,c=c,texture=texture,mass=mass,label=label,draw=draw,move=move,body_list=body_list)
        self.order=4
    def update_vel(self,params_vel_pos=None,dt=0): #kutta4 velocities
        if self.move:
            self.v0=self.v
            k1=params_vel_pos[0]['copy_'+self.label]['a']
            k2=params_vel_pos[1]['copy_'+self.label]['a']
            k3=params_vel_pos[2]['copy_'+self.label]['a']
            k4=params_vel_pos[3]['copy_'+self.label]['a']
            self.v += (k1 + k2*2 + k3*2 + k4)*dt/6
