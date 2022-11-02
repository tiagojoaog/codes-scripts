#!/usr/bin/env python3

from vpython import *
from matplotlib import pyplot as plt
import numpy as np
import sys,os,imageio,time,pickle,gc
from os.path import expanduser

try:
    import tty, termios, fcntl
    fd = sys.stdin.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
except:
    import keyboard

__all__ = ['orbit','center_of_mass','hour_min_sec','getch']


def center_of_mass(object_list):  #outputs the position of the center of mass, along with the collective mass, respectively.
    M=0
    sum_xm=vector(0,0,0)
    for i in object_list:
        x=i.object.pos
        m=i.mass
        sum_xm+=x*m
        M+=m
    return sum_xm/M,M

def hour_min_sec(time):   # converts total seconds to hh:mm:ss format
    h=time//3600
    m=(time % 3600)//60
    s=((time %3600)%60)//1
    return str(int(h))+':'+str(int(m))+':'+str(int(s))

def getch():                                    # this definition check what key has been pressed. It is only called on Linux and Mac OS. 
    fd = sys.stdin.fileno()                     # for windows, keyboard is used instead.
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        time.sleep(2E-4)
        ch = sys.stdin.read(1)
    except:
        ch = None
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch


class orbit(): # class responsible for iterating over bodies. Has all the tools required
    def __init__(self, G= 2.959122083e-4,Ke=0,draw=True,frames=False,traj_show_save='show',dt=4.166666e-2,RCM=vector(1,0,0),capture=False,is_adaptive=True,bodies=[]):
        if draw: #animates
            from vpython.no_notebook import stop_server
            global stop_server
            scene.caption = """Planetary Motion:
            To rotate "camera", drag with right button or Ctrl-drag.
            To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
                On a two-button mouse, middle is left + right.
            To pan left/right and up/down, Shift-drag.
            Touch screen: pinch/extend to zoom, swipe or two-finger rotate.
            Hold q to quit."""
        
            scene.forward = vector(0,-.3,-1)
        print("Planetary motion initiated.")
        print("-- Hold q to terminate\n")
        self.G              = G
        self.Ke             = Ke
        self.frames         = frames
        self.dt             = dt
        self.capture        = capture
        self.traj_show_save = traj_show_save 
        self.is_adaptive    = is_adaptive
        self.bodies         = bodies
        self.draw           = draw
        if self.bodies:
            self.order      = bodies[0].order
        self.RCM            = RCM
        self.t              = 0    
        self.n              = 0               
        self.filenames      = []
    def finish(self): #Final print
        print('Years:',"{:.2f}".format(self.t/(365.25)),' dt(days):',"{:.6f}".format(self.dt))
        print('Terminated by pressing q.')
        if self.capture:
            print('Waiting for PNG to be stored...')
            sleep(3)
    def close(self): #closes the program
        if self.draw:
            stop_server()
        sys.exit(1)
    def flush_cache(self,save,string): #output properties, energies, etc... to a pickle file
        with open(string+'.pkl', 'ab') as pickel:
            pickle.dump(save,pickel)
    def read_cached_quantities(self,string): #reads pickle files for energies, angular momentum, step (dt)
        cached = None
        if os.path.exists(string+'.pkl'):
            with open(string+'.pkl', 'rb') as pickel:
                in_range=True
                while in_range:
                    try:
                        loaded_pickle = pickle.load(pickel)
                        for i in loaded_pickle:
                            if not cached:
                                cached = [i]
                            else:
                                cached.append(i)
                    except:
                        in_range=False
        return cached
    def read_cached_properties(self,string): #reads the object's properties from pickle files
        cached = None
        if os.path.exists(string+'.pkl'):
            with open(string+'.pkl', 'rb') as pickel:
                in_range=True
                iterate=0
                while in_range:
                    try:
                        loaded_pickle = pickle.load(pickel)
                        for i in loaded_pickle:
                            iterate2=0
                            for j in loaded_pickle[i]:
                                if not cached:
                                    cached = {}
                                    cached[i] = {}
                                    cached[i]['pos']   = [j['pos']]
                                    cached[i]['v']     = [j['v']]
                                    
                                else:
                                    if iterate==0 and iterate2==0:
                                        cached[i] = {}
                                        cached[i]['pos']   = [j['pos']]
                                        cached[i]['v']     = [j['v']]
                                    else:
                                        cached[i]['pos']   += [j['pos']]
                                        cached[i]['v']     += [j['v']]
                                    
                                cached[i]['q']     = j['q']    
                                cached[i]['mass']  = j['mass'] 
                                cached[i]['draw']  = j['draw'] 
                                cached[i]['move']  = j['move']
                                cached[i]['light']  = j['light']
                                cached[i]['color']  = j['color']
                                cached[i]['texture']  = j['texture']
                                cached[i]['radius']  = j['radius']
                                iterate2+=1
                        iterate+=1
                    except:
                        in_range=False
        return cached
    def clear_cache(self): #deletes pickle files
        here = os.getcwd()
        files = os.listdir(here)
        for i in files:
            if i.split('.')[-1]=='pkl':
                os.remove(os.path.join(here,i))
    def update_cache(self,string,parameter): #updates pickle files
        self.flush_cache(parameter,string)
    def iterate_cache(self,KE=None,PE=None,totE=None,step=None,ang_mom=None): #uses definitions created above to update all properties 
        self.update_cache('kinetic_energy',KE)
        self.update_cache('potential_energy',PE)
        self.update_cache('total_energy',totE)
        self.update_cache('angular_momentum',ang_mom)
        self.update_cache('steps',step)
        properties={}
        for i in self.bodies:
            properties[i.label]=i.images
            i.clear()
        self.update_cache('properties',properties)
        del properties
    def iterate(self,bodies=[],RCM=vector(0,0,0),params_vel_pos=None,acel=None,dt=0,order=1,actual_run=True): #calculates accelerations,velocities,positions and kutta coefficients
        vel_pos_calc={}                                                                                       
        for i in bodies:
            vel_pos_calc[i.label]={}
            if not acel:
                i.a=vector(0,0,0)
            for j in bodies:
                if i!=j:
                    if not acel:
                        i.acceleration(b_object=j,G=self.G,Ke=self.Ke)
            if not actual_run or order==1:
                i.update_vel_euler(dt=dt)
            else:
                i.update_vel(params_vel_pos=params_vel_pos,dt=dt)
            i.update_pos(dt=dt)
            vel_pos_calc[i.label]['pos']=vector(i.object.pos)
            vel_pos_calc[i.label]['v']=vector(i.v)
            i.a=vector(0,0,0)
            for j in bodies:
                if i!=j:
                    i.acceleration(b_object=j,G=self.G,Ke=self.Ke)
            vel_pos_calc[i.label]['a']=i.a                  # Uses updated acceleration to calculate kn on runge-kutta methods
        KE=0
        PE=0
        ANG_MOM=vector(0,0,0)
        for i in self.bodies:
            i.energies(self.bodies,G=self.G,Ke=self.Ke,RCM=RCM)
            KE+=i.KE
            PE+=i.PE
            ANG_MOM+=i.ANG_MOM
        ANG_MOM=mag(ANG_MOM)
        PE=PE/2                #double-counting fix
        return KE,PE,ANG_MOM,vel_pos_calc
    def fake_run(self,body_list=None,order=1,iterations=1,dt=0): #uses iterate to perform a fake run (no actual positions are updated)
        copied_bodies=[]                                         #this is used for adaptive step calculation and kutta coefficients
        etot_list=[]
        ANG_mom_list=[]
        iter_vel_pos_calc=[]
        if self.order==4:
            dt=dt/2
            iterations=3

        for i in body_list:
            copied_object=i.copy(label='copy_'+i.label)
            copied_bodies.append(copied_object)
        vel_pos_calc={}
        for i in copied_bodies:
            i.a=vector(0,0,0)
            vel_pos_calc[i.label]={}
            for j in copied_bodies:
                if i!=j:
                    i.acceleration(b_object=j,G=self.G,Ke=self.Ke)
                    vel_pos_calc[i.label]['pos']=i.object.pos
                    vel_pos_calc[i.label]['v']=i.v
                    vel_pos_calc[i.label]['a']=i.a
        iter_vel_pos_calc.append(vel_pos_calc)              
        for cyc in range(iterations):
            RCM,dummy_mass=center_of_mass(copied_bodies)
            del dummy_mass
            if cyc == 1 and order == 4:
                for i in copied_bodies:
                    i.object.pos = iter_vel_pos_calc[0][i.label]['pos']
                    i.v = iter_vel_pos_calc[0][i.label]['v']
                    i.a = iter_vel_pos_calc[1][i.label]['a']
                KE,PE,ANG_MOM,vel_pos_calc=self.iterate(bodies=copied_bodies,RCM=RCM,order=1,acel=True,dt=dt,actual_run=False)
            elif cyc == 2 and order == 4:
                for i in copied_bodies:
                    i.object.pos = iter_vel_pos_calc[0][i.label]['pos']
                    i.v = iter_vel_pos_calc[0][i.label]['v']
                    i.a = iter_vel_pos_calc[2][i.label]['a']
                KE,PE,ANG_MOM,vel_pos_calc=self.iterate(bodies=copied_bodies,RCM=RCM,order=1,acel=True,dt=2*dt,actual_run=False)
            else:
                KE,PE,ANG_MOM,vel_pos_calc=self.iterate(bodies=copied_bodies,RCM=RCM,order=1,dt=dt,actual_run=False)
            Etot=KE+PE
            etot_list.append(Etot)
            ANG_mom_list.append(ANG_MOM)
            iter_vel_pos_calc.append(vel_pos_calc)
        del copied_bodies,copied_object,vel_pos_calc,KE,PE,ANG_MOM
        return etot_list,ANG_mom_list,iter_vel_pos_calc
    def adaptive(self,body_list,threshold_adaptive=0,max_dt=0): #adaptive step based on the conservation of angular momentum and total energy
        self.dt=max_dt                                          #if error below threshold, dt is used for the actual run.
        while True:
            etot_list,ANG_mom_list,dummy_dict=self.fake_run(body_list=body_list,dt=self.dt,iterations=2)
            del dummy_dict
            error=[]
            old_num1=None
            for num1,num2 in zip(etot_list,ANG_mom_list):
                if not old_num1:
                    old_num1=num1
                    old_num2=num2
                else:
                    if num1==0 and old_num1==0:
                        error.append(0)
                    elif old_num1==0:
                        error.append(abs((num1-old_num1)/num1))
                    else:
                        error.append(abs((num1-old_num1)/old_num1))
                    
                    if num2==0 and old_num2==0:
                        error.append(0)
                    elif old_num2==0:
                        error.append(abs((num2-old_num2)/num2))
                    else:
                        error.append(abs((num2-old_num2)/old_num2))
                    old_num1=num1
                    old_num2=num2
            if max(error)<threshold_adaptive:
                del error,old_num1,old_num2,num1,num2,etot_list,ANG_mom_list
                return True
            else:
                self.dt=self.dt/1.5
                if self.dt < 5E-5:
                    print("\n\nTimestep too small:",self.dt)
                    print("  Your objects most likely collided.") 
                    return False
    def make_gif(self,capture,filenames): #captures a gif from animation to use in presentations
        if capture:
            home = expanduser("~")
            path=home+'/Downloads'
            here=os.getcwd()

            for i in filenames:
                os.system('mv '+path+'/'+i+'.png .')
            print('making GIF...')
            with imageio.get_writer(here+'/movie.gif', mode='I') as writer:
                for filename in filenames:
                    filename+='.png'
                    image = imageio.imread(filename)
                    writer.append_data(image)
            os.system('rm image_*png')
            print('done.')

    def run_animation(self,properties,timeperiod):     #This is used in case visualization of the animation is preferred after the calculation. In this case, the animation is not
        from vpython.no_notebook import stop_server    #influenced by the calculation. Animation is done based on the information of the pickle files
        scene.caption = """Planetary Motion:
        To rotate "camera", drag with right button or Ctrl-drag.
        To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
            On a two-button mouse, middle is left + right.
        To pan left/right and up/down, Shift-drag.
        Touch screen: pinch/extend to zoom, swipe or two-finger rotate.
        Press q to quit."""

        scene.forward = vector(0,-.3,-1)
        
        for i in properties:
            c       = properties[i]['color']
            r       = properties[i]['radius']
            texture = properties[i]['texture']
            properties[i]['object']   = sphere(pos=properties[i]['pos'][0], radius=r*4.65047e-3, color=c, texture=texture,
                                        make_trail=True, trail_type='points', interval=10, retain=50)
            if properties[i]['light']:
                properties[i]['local_light'] = local_light(pos=properties[i]['pos'][0],
                                               color=c)
        n=0
        if self.capture:
            name='image_'+str(self.t)
            scene.capture(name)
            self.filenames.append(name)
        for i in timeperiod:
            if self.capture and not self.frames:
                print('Cannot capture images if frames are not stated.')
                break
            if self.frames:
                rate(self.frames)
                if (n%(self.frames/4))==0:
                    name='image_'+str(i)
                    if self.capture:
                        scene.capture(name)
                    self.filenames.append(name)
            for j in properties:
                if properties[j]['move']:
                    properties[j]['object'].pos = properties[j]['pos'][n]
                    if properties[j]['light']:
                        properties[j]['local_light'].pos = properties[j]['object'].pos
            n+=1
            print('Years:',"{:.2f}".format(i/(365.25)),end='\r')
            if 'termios' in sys.modules:
                char = getch()
            else:
                if keyboard.is_pressed('q'):
                    char='q'
                else:
                    char=None

            if (char == 'q'):
                self.finish()
                break
        self.make_gif(self.capture,self.filenames)
    def analysis(self,traj_show_save=None,animate=False):                  #creates energy, angular momentum and position graphs. Animation definition is called here.
            self.step_list          = self.read_cached_quantities('steps')
            self.KE_list            = self.read_cached_quantities('kinetic_energy')
            self.PE_list            = self.read_cached_quantities('potential_energy')
            self.totE_list          = self.read_cached_quantities('total_energy')
            self.ang_mom_list       = self.read_cached_quantities('angular_momentum')
            self.properties         = self.read_cached_properties('properties')
            self.step_list.pop(0)
            plt.plot(self.step_list,self.KE_list,label='Kinetic')
            plt.plot(self.step_list,self.PE_list,label='Potential')
            plt.plot(self.step_list,self.totE_list,label='Total')
            plt.legend(loc='best')
            plt.savefig('energies.pdf')
            plt.close()
            plt.plot(self.step_list,self.totE_list,label='Total')
            plt.savefig('mechanical_energy.pdf')
            plt.close()
            plt.plot(self.step_list,self.ang_mom_list,label='Angular Momentum')
            plt.legend(loc='best')
            plt.savefig('angular_momentum.pdf')
            plt.close()
            fig=plt.figure()
            ax=plt.axes(projection='3d')
            for i in self.properties:
                positions = self.properties[i]['pos']
                X=[]
                Y=[]
                Z=[]
                for j in range(len(positions)):
                    X.append(self.properties[i]['pos'][j].x)
                    Y.append(self.properties[i]['pos'][j].y)
                    Z.append(self.properties[i]['pos'][j].z)
                ax.plot3D(X,Y,Z,label=i)
            ax.view_init(155,-95)
            plt.legend(loc='best')
            if animate:
                self.run_animation(self.properties,self.step_list)  
            if traj_show_save == 'show':
                plt.show()
            else:
                plt.savefig('trajectories.pdf')
                plt.close()
    def run(self,threshold_adaptive=1e-4,max_dt=1,timelimit=None,clear_cache=False): #Main definition. Where the magic happens. Adaptive step calculated, coefficients calculated,
        start = time.time()                                                          #properties updated (energies, pos, velocities, etc...)
        self.threshold_adaptive=threshold_adaptive
        self.max_dt=max_dt
        KE_list=[]
        PE_list=[]
        totE_list=[]
        step_list=[0]
        ang_mom_list=[]
        if self.capture:
            name='image_'+str(self.t)
            scene.capture(name)
            self.filenames.append(name)
        while True:
            if self.capture and not self.frames:
                print('Cannot capture images if frames are not stated.')
                break
            if self.frames:
                rate(self.frames)
                if (self.n%(self.frames/4))==0:
                    name='image_'+str(self.t)
                    if self.capture:
                        scene.capture(name)
                    self.filenames.append(name)
            if self.is_adaptive:
                continuation=self.adaptive(self.bodies,threshold_adaptive=self.threshold_adaptive,max_dt=self.max_dt)
                if not continuation:
                    break
            dummy1,dummy2,params_vel_pos=self.fake_run(body_list=self.bodies,dt=self.dt,order=self.order)
            KE,PE,ANG_MOM,dummy3=self.iterate(bodies=self.bodies,RCM=self.RCM,params_vel_pos=params_vel_pos,order=self.order,dt=self.dt,actual_run=True)
            self.t+=self.dt
            self.n+=1
            totE=KE+PE
            KE_list.append(KE)
            PE_list.append(PE)  
            totE_list.append(totE)
            ang_mom_list.append(ANG_MOM)
            step_list.append(self.t)
            if (self.n%500)==0:           #Dump cache and clear memory after each 500 iterations
                self.iterate_cache(KE=KE_list,PE=PE_list,totE=totE_list,step=step_list,ang_mom=ang_mom_list)
                KE_list=[]
                PE_list=[]
                totE_list=[]
                step_list=[]
                ang_mom_list=[]
                gc.collect()
            self.RCM,dummy_mass=center_of_mass(self.bodies)
            del dummy_mass,dummy1,dummy2,dummy3
            del params_vel_pos,KE,PE,totE,ANG_MOM
            print('Years:',"{:.2f}".format(self.t/(365.25)),' dt(days):',"{:.6f}".format(self.dt),end='\r')
            if 'termios' in sys.modules:
                char = getch()
            else:
                if keyboard.is_pressed('q'):
                    char='q'
                else:
                    char=None
            if (char == 'q'):
                self.finish()
                break
            if timelimit:
                if (self.t)/365.25 > timelimit:
                    self.finish()
                    break
        if step_list:
            self.iterate_cache(KE=KE_list,PE=PE_list,totE=totE_list,step=step_list,ang_mom=ang_mom_list)
        end=time.time()
        print("\n Runtime:",hour_min_sec(end-start))
        self.make_gif(self.capture,self.filenames)
        self.analysis(traj_show_save=self.traj_show_save)
        if clear_cache:
            self.clear_cache()
