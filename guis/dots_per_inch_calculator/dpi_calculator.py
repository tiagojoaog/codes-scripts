#!/usr/bin/env python3

from functools import partial
from tkinter import messagebox

import tkinter as tk
import tkinter.ttk as ttk
import numpy as np
#import os,sys

def converter(which='cm_to_inches',value=0):
    if which == 'cm_to_inches':
      ratio = 0.393701
    
    elif which == 'inches_to_cm':
      ratio = 1/0.393701

    elif which == 'inches_to_inches':
      ratio = 1
    
    else:
      raise NameError(f"Conversion '{which}' not recognized.")

    try:
      return float(value)*ratio

    except ValueError:
      if not value:
         value=None
      else:
         return value


class DPI_calculator():
    def __init__(self):
        self.windows = {}
        self.dpis    = []
        self.Window(x=590,y=460,title='DPI Calculator GUI',window_name='root')
        self.Label(window_name='root',label='Resolution', text='Window Resolution\n')
        self.Entry(window_name='root',label='Width',text='Width',width=10)
        self.Entry(window_name='root',label='Height',text='Height',width=10)
        self.Button(window_name='root',label='Open_topwindow',text='Open',width=19, command=partial(self.Window, title='Calculate DPI', window_name='DPI'))
        self.Button(window_name='root',label='Close_all',text='Close',width=5, command=partial(self.destroy))

        self.Label(window_name='root',label='instructions', text=''' 

    INSTRUCTIONS:

    To evaluate your mouse's DPI, do a control+click(LMB) on either side of the white screen
    and another on the other end. Have a ruler with you to measure how much your mouse 
    (hardware) has moved on your pad. Then insert how much it moved on the 
    'Mouse Displacement' box and click 'Calculate'. 

    Repeat this as many times as you want, the numbers will be stored.
    
    To obtain the average of all the values obtained, click on the 'Average' button.
    This will provide the average, maximum value, minimum value and the standard
    deviation.
    
    'Clear DPI values' will wipe all the numbers obtained to reset the evaluation.
    
    Window Resolution is not relevant, it is just a set of parameters to define 
    how big you want your window to be.

    HINT: Make your window as large as possible so that you will have more space
              to move your mouse, leading to a more accurate DPI calculation.

'''


)
          
    def Window(self, *, x=500, y=500, title, window_name, event=None):
        if window_name in self.windows:
            self.destroy(window=window_name)

        self.windows[window_name] = {}

        if window_name == 'root':
            self.windows['root']['window'] = tk.Tk()
            self.windows[window_name]['window'].geometry(f'{x}x{y}')
            self.windows[window_name]['window'].title(title)

        else:
            self.windows[window_name]['window'] = tk.Toplevel()
            self.first_value = None
            self.which       = None
            self.size_window_x=self.windows['root']['entry']['Width'].get()
            self.size_window_y=self.windows['root']['entry']['Height'].get()
            self.windows[window_name]['window'].geometry(f'{self.size_window_x}x{self.size_window_y}')
            self.windows[window_name]['window'].title(title)
            self.windows[window_name]['canvas'] = tk.Canvas(self.windows[window_name]['window'], bg="white", height=int(self.size_window_y)-220, width=int(self.size_window_x)-30, relief='sunken', borderwidth=1)
            self.windows[window_name]['canvas'].grid(row=1,column=1,padx=10,columnspan=3)
            self.windows[window_name]['mouse'] = {}
            self.windows[window_name]['mouse']['mouse_coordinates'] = ttk.Label(self.windows[window_name]['window'],text=f' Width={0:>6} Height={0:>6} ',background="white", relief='sunken', borderwidth=2)
            self.windows[window_name]['mouse']['refresh'] = 0
            self.windows[window_name]['canvas'].bind('<Motion>',partial(self.mouse_coordinates,window_name))
            self.windows[window_name]['canvas'].bind('<Control-Button-1>',partial(self.save_value,window_name))
            self.windows[window_name]['mouse']['mouse_coordinates'].grid(row=20,column=1,columnspan=10)
            self.Entry(window_name=window_name,width=10,label='mouse_movement',text='Mouse Displacement')
            self.windows[window_name]['label']['mouse_movement'].grid(row=2,column=1,columnspan=2)
            self.windows[window_name]['entry']['mouse_movement'].grid(row=3,column=1,sticky=tk.E)
            self.windows[window_name]['options'] = {}
            self.windows[window_name]['options']['cm_inches_var']  = tk.StringVar(self.windows[window_name]['window'])
            self.windows[window_name]['options']['cm_inches_drop'] = tk.OptionMenu(self.windows[window_name]['window'], self.windows[window_name]['options']['cm_inches_var'], 'cm', 'Inches',command=partial(self.which_one,window_name))
            self.windows[window_name]['options']['cm_inches_var'].set('cm')
            self.which_one(window_name)
            self.windows[window_name]['options']['cm_inches_drop'].grid(row=3,column=2,sticky=tk.W)
            self.Button(window_name=window_name,label='calculate',text='Calculate',width=20,command=partial(self.calculate,window_name))
            self.windows[window_name]['button']['calculate'].grid(row=5,column=1,columnspan=2)
            self.windows[window_name]['final_result'] = ttk.Label(self.windows[window_name]['window'],text=f' DPI={0:>5}, {len(self.dpis):>3} number(s) saved ', background="white", relief='sunken', borderwidth=2)
            self.windows[window_name]['final_result'].grid(row=5,column=3)
            self.Button(window_name=window_name,label='clear',text='Clear DPI values (reset average)',width=25,command=partial(self.clear, window=window_name,clear_all=True))
            self.Button(window_name=window_name,label='statistics',text='Average', width=20,command=partial(self.statistical_analysis, window_name))
            self.windows[window_name]['button']['clear'].grid(row=6,column=3)
            self.windows[window_name]['button']['statistics'].grid(row=6,column=1,columnspan=2)
            self.windows[window_name]['whitespace1'] = ttk.Label(self.windows[window_name]['window'],text=' ').grid(row=8,column=1)
            self.windows[window_name]['average_dpi'] = ttk.Label(self.windows[window_name]['window'],text=f' Average DPI:{0:>5.3f} with a standard deviation of {0:>5.3f},\n maximum value of {0:>5.3f} and minimum value of {0:>5.3f} ', background="white", relief='sunken', borderwidth=2)
            self.windows[window_name]['average_dpi'].grid(row=9,column=1,columnspan=5)
            self.windows[window_name]['whitespace1'] = ttk.Label(self.windows[window_name]['window'],text=' ').grid(row=11,column=1)

    def Widgets(self, *, window_name, widget_type, widget, label):
        if widget_type not in self.windows[window_name]:
            self.windows[window_name][widget_type] = {}

        self.windows[window_name][widget_type][label] = widget

    def Label(self, *, window_name, label, text):
        self.Widgets(widget_type='label', window_name=window_name, widget=ttk.Label(self.windows[window_name]['window'],text=text), label=label)

    def Entry(self, justify=tk.CENTER, *, window_name, label, text, width):
        self.Label(window_name=window_name, label=label, text=text)
        self.Widgets(widget_type='entry', widget=ttk.Entry(self.windows[window_name]['window'],width=width,justify=justify), window_name=window_name, label=label)
   
    def Button(self, *, window_name, label, text, width, command):
        self.Widgets(widget_type='button', widget=ttk.Button(self.windows[window_name]['window'],width=width, text=text, command=command), window_name=window_name, label=label)

    def clear(self, window=None ,clear_all=False, event=None):
        self.first_value  = None
        self.second_value = None

        self.windows[window]['canvas'].destroy()
        self.windows[window]['canvas'] = tk.Canvas(self.windows[window]['window'], bg="white", height=int(self.size_window_y)-220, width=int(self.size_window_x)-30, relief='sunken', borderwidth=1)
        self.windows[window]['canvas'].grid(row=1,column=1,padx=10,columnspan=3)
        self.windows[window]['canvas'].bind('<Motion>',partial(self.mouse_coordinates,window))
        self.windows[window]['canvas'].bind('<Control-Button-1>',partial(self.save_value,window))

        if clear_all:
            self.dpis = []
            self.windows[window]['average_dpi'].destroy()
            self.windows[window]['final_result'].destroy()
            self.windows[window]['average_dpi'] = ttk.Label(self.windows[window]['window'],text=f' Average DPI:{0:>5.3f} with a standard deviation of {0:>5.3f},\n maximum value of {0:>5.3f} and minimum value of {0:>5.3f} ', background="white", relief='sunken', borderwidth=2)
            self.windows[window]['average_dpi'].grid(row=9,column=1,columnspan=5)
            self.windows[window]['final_result'] = ttk.Label(self.windows[window]['window'],text=f' DPI={0:>5d}, {len(self.dpis):>3} number(s) saved ', background="white", relief='sunken', borderwidth=2)
            self.windows[window]['final_result'].grid(row=5,column=3)

    def statistical_analysis(self, window_name):
        mean    = np.mean(self.dpis)
        std     = np.std(self.dpis)
        max_val = max(self.dpis) 
        min_val = min(self.dpis)
        self.windows[window_name]['average_dpi'] = ttk.Label(self.windows[window_name]['window'],text=f' Average DPI:{mean:>5.3f} with a standard deviation of {std:>5.3f},\n maximum value of {max_val:>5.3f} and minimum value of {min_val:>5.3f} ', background="white", relief='sunken', borderwidth=2)
        self.windows[window_name]['average_dpi'].grid(row=9,column=1,columnspan=5)
        print(f'Mean DPI is {mean} with a standard deviation of {std}.')

    def calculate(self,window_name,event=None):
        inches = converter(self.which, self.windows[window_name]['entry']['mouse_movement'].get())

        if not isinstance(inches,float):
            messagebox.showerror('ERROR!',f"'{inches}' not recognized.\nPlease use a valid number to define how much your physical mouse has moved either in cm or in inches.")
        else:
            dpi = abs(round((int(self.second_value)-int(self.first_value))/inches,3))
            self.dpis.append(dpi)
            self.clear(window=window_name)
            self.windows[window_name]['final_result'].destroy()
            self.windows[window_name]['final_result'] = ttk.Label(self.windows[window_name]['window'],text=f' DPI={dpi:>5.3f}, {len(self.dpis):>3} number(s) saved ', background="white", relief='sunken', borderwidth=2)
            self.windows[window_name]['final_result'].grid(row=5,column=3)
            print(f' option = {self.which} inches={inches}   width displacement (pixels) = {dpi*inches}')
            print(' DPI List:'+self.list_to_string(self.dpis))

    def mouse_coordinates(self, window, event):
        x, y = event.x, event.y
       
        self.windows[window]['mouse']['refresh'] +=1
        self.windows[window]['mouse']['mouse_coordinates'] = ttk.Label(self.windows[window]['window'],text=f' Width={x:>6} Height={y:>6} ', background='white', relief='sunken', borderwidth=2)
        self.coordinates = [x,y]
        self.windows[window]['mouse']['mouse_coordinates'].grid(row=20,column=1,padx=0,columnspan=10)

    def save_value(self,window_name,event=None):
        self.make_circle(self.coordinates[0],self.coordinates[1],3,self.windows[window_name]['canvas'])
        if not self.first_value:
            canvas_label = 'coord1'
            self.first_value = self.coordinates[0]
        else:
            canvas_label = 'coord2'
            self.second_value = self.coordinates[0]
        ttk.Label(self.windows[window_name]['canvas'],text=f' x={self.coordinates[0]:>4}\n y={self.coordinates[1]:>4} ',background='white').place(x=self.coordinates[0]-25,y=self.coordinates[1]+5)
        print(f'coordinate x = {self.coordinates[0]}')

    def make_circle(self,x,y,r,canvas):
        x0=x-r
        y0=y-r
        x1=x+r
        y1=y+r
        return canvas.create_oval(x0,y0,x1,y1,fill='red')

    def which_one(self,window_name,event=None):
        variable = self.windows[window_name]['options']['cm_inches_var'].get()             
        if variable == 'cm':
           self.which = 'cm_to_inches'

        else:
           self.which = 'inches_to_inches'

    def list_to_string(self,_list):
        _string = ''
        for i in _list:
            _string+=f' {i}'
        return _string

    def destroy(self,window=None, event=None):
        if window:
            self.windows[window]['window'].destroy()
            self.windows[window] = None  
        else:
            for window in self.windows:
                try:
                    self.windows[window]['window'].destroy()
                except tk.TclError:
                    pass
                self.windows[window] = None             

    def Loop(self):
        self.windows['root']['label']['Resolution'].grid(row=0,column=2,columnspan=3)
        self.windows['root']['label']['Width'].grid(row=1,column=3)
        self.windows['root']['label']['Height'].grid(row=1,column=4,sticky=tk.W)
        self.windows['root']['entry']['Width'].grid(row=2,column=3,sticky=tk.E)
        self.windows['root']['entry']['Height'].grid(row=2,column=4,sticky=tk.W)
        self.windows['root']['button']['Open_topwindow'].grid(row=3,column=3,columnspan=2)
        self.windows['root']['label']['instructions'].grid(row=4,column=0,columnspan=10)
        self.windows['root']['button']['Close_all'].grid(row=5,column=3,columnspan=2)
        self.windows['root']['window'].mainloop()


if __name__ == '__main__' :
    DPI = DPI_calculator()
    DPI.Loop()
