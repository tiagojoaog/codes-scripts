#!/usr/bin/env python3

import pickle
import numpy as np
import tkinter as tk
#import tkinter.ttk as ttk
from tkinter import messagebox
from functools import partial
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from fancy_plots import fancy_plots
from PIL import ImageTk, Image

import os,sys

'''Fancy plots GUI'''

windows   = {}
execution = {}

add_text=[]
parameters=[]
dict_list=[]   #combined paths



''' Tooltips '''

class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                      background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

''''''

''' Cut, Copy and Paste: Menu'''

def make_textmenu(root):
	the_menu = tk.Menu(root, tearoff=0)
	the_menu.add_command(label="Cut")
	the_menu.add_command(label="Copy")
	the_menu.add_command(label="Paste")
	the_menu.add_separator()
	the_menu.add_command(label="Select all")
	return the_menu
def callback_select_all(root,event):
	# select text after 50ms
	root.after(50, lambda:event.widget.select_range(0, 'end'))

def show_textmenu(the_menu,event):
	e_widget = event.widget
	the_menu.entryconfigure("Cut",command=lambda: e_widget.event_generate("<<Cut>>"))
	the_menu.entryconfigure("Copy",command=lambda: e_widget.event_generate("<<Copy>>"))
	the_menu.entryconfigure("Paste",command=lambda: e_widget.event_generate("<<Paste>>"))
	the_menu.entryconfigure("Select all",command=lambda: e_widget.select_range(0, 'end'))
	the_menu.tk.call("tk_popup", the_menu, event.x_root, event.y_root)


''''''


def add_path(entry,dict_list):
    dict_list.append(entry)

def mouse_coordinates(root,execution,ax1,param_dict,event):
    x, y = event.x, event.y
    if isinstance(param_dict,dict):
      xfactor = float(param_dict['boxsize'].get().split(',')[0])
      yfactor = float(param_dict['boxsize'].get().split(',')[1])
    elif isinstance(param_dict,list):
      for item in param_dict:
         if item.split('=')[0]=='boxsize':
           xfactor = float(item.split('=')[1].split(',')[0])
           yfactor = float(item.split('=')[1].split(',')[1])
    cursor_xmin=15*xfactor
    cursor_xmax=108*xfactor
    cursor_ymin=107*yfactor
    cursor_ymax=15*yfactor
    if int(x) < cursor_xmin:
        x=cursor_xmin
    elif int(x) > cursor_xmax:
        x=cursor_xmax
    if int(y) < cursor_ymax:
        y=cursor_ymax
    elif int(y) >cursor_ymin:
        y=cursor_ymin


    xmin,xmax = ax1.get_xlim()
    ymin,ymax = ax1.get_ylim()

    x = xmin+(x-cursor_xmin)*(xmax-xmin)/(cursor_xmax-cursor_xmin)
    y = ymin+(y-cursor_ymin)*(ymax-ymin)/(cursor_ymax-cursor_ymin)
    execution['initial_coordinates'].destroy()
    execution['coordinates'] = tk.Label(root,text=f'x={x:.3f} y={y:.3f}')
    execution['coordinates'].grid(row=10,column=1,sticky=tk.S)


def convert_path_to_list(paths_dict,full_pathway):
    full_list      = []
    paths          = []
    adict          = {}
    labelling      = []

    for name in paths_dict:
        name = name.split('_')[0].split('-')[0]
        if name not in paths:
            paths.append(name)
    paths.sort()
    for path in paths:
        for name in full_pathway:
            full_name = f'{path}_{name}'
            if full_name in paths_dict:
                adict[name] = paths_dict[full_name]
        full_list.append(adict)
        adict = {}
        legendlabel = f'{path}-legendlabel'
        if legendlabel in paths_dict:
            labelling.append(paths_dict[legendlabel])
    return full_list,labelling

def update_graph(param_dict,parameters,paths_dict,full_pathway,root,matplotlib_params,add_text_dict,add_text):
    matplotlib_params['ax1'],matplotlib_params['fig']=generate_figure(param_dict,parameters,paths_dict,full_pathway,add_text_dict,add_text)
    if 'canvas' in matplotlib_params:
        matplotlib_params['canvas']._tkcanvas.destroy()
    matplotlib_params['canvas'] = FigureCanvasTkAgg(matplotlib_params['fig'],master=root)
    matplotlib_params['canvas'].draw()
    matplotlib_params['canvas'].get_tk_widget().grid(row=9,column=1,padx=30)
    if not matplotlib_params['ax1'] and not matplotlib_params['fig']:
        return True

def obtain_boxsizes(parameters):
    for word in parameters:
           if word.split('=')[0] == 'boxsize':
              scale_window_x = float(word.split('=')[1].split(',')[0])
              scale_window_y = float(word.split('=')[1].split(',')[1])
              if scale_window_x <2:
                 scale_window_x = 2
    return scale_window_x,scale_window_y

def tight_layout_on_press(param_dict,parameters,paths_dict,full_pathway,matplotlib_params,add_text_dict,add_text,windows,event=None):
    if windows['tight_layout']:
        quit_window(windows,'tight_layout',matplotlib_params,'canvas_tl')
    windows['tight_layout'] = tk.Toplevel()
    scale_window_x,scale_window_y=obtain_boxsizes(parameters)
    windows['tight_layout'].geometry(f'{int((120*scale_window_x)+80)}x{int((120*scale_window_y)+80)}')
    windows['tight_layout'].title('Fancy Plots - Preview')
    fig = matplotlib_params['fig']
    fig.tight_layout()
    if 'canvas_tl' in matplotlib_params:
        matplotlib_params['canvas_tl']._tkcanvas.destroy()
    matplotlib_params['canvas_tl'] = FigureCanvasTkAgg(fig,master=windows['tight_layout'])
    matplotlib_params['canvas_tl'].draw()
    matplotlib_params['canvas_tl'].get_tk_widget().grid(row=9,column=1,padx=30)
    windows['tight_layout'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'tight_layout',matplotlib_params,'canvas_tl'))

def tight_layout_on_release(windows,matplotlib_params,event=None):
    quit_window(windows,'tight_layout',matplotlib_params,'canvas_tl')

def graph(param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text):
    if windows['graph_window']:
        try:
            destroy=update_graph(param_dict,parameters,paths_dict,full_pathway,windows['graph_window'],matplotlib_params,add_text_dict,add_text)
            execution['initial_coordinates'] = tk.Label(windows['graph_window'],text=f'x=0.000 y=0.000')
            execution['initial_coordinates'].grid(row=10,column=1,sticky=tk.S)
            windows['graph_window'].bind('<Control-Motion>',partial(mouse_coordinates,windows['graph_window'],execution,matplotlib_params['ax1'],parameters))
        except tk.TclError:
            windows['graph_window'] = None
            destroy = graph(param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text)
        if destroy:
            quit_window(windows,'graph_window',matplotlib_params,'canvas')
        else:
            scale_window_x,scale_window_y=obtain_boxsizes(parameters)
            windows['graph_window'].geometry(f'{int((120*scale_window_x)+80)}x{int((120*scale_window_y)+80)}')
    else:
        windows['graph_window'] = tk.Toplevel()
        make_parameters(param_dict,parameters)
        scale_window_x,scale_window_y=obtain_boxsizes(parameters)
        windows['graph_window'].geometry(f'{int((120*scale_window_x)+80)}x{int((120*scale_window_y)+80)}')
        windows['graph_window'].title('Fancy Plots - Graph')
        
        make_parameters(param_dict,parameters)

        execution = {}
        execution['initial_coordinates'] = tk.Label(windows['graph_window'],text=f'x=0.000 y=0.000')
        execution['initial_coordinates'].grid(row=10,column=1,sticky=tk.S)
        
        cleanup_button = tk.Button(windows['graph_window'],text='Hold to Preview Graph with Tight Layout')
        cleanup_button.grid(row=11,column=1,sticky=tk.S)
        cleanup_button.bind('<ButtonPress-1>',partial(tight_layout_on_press,param_dict,parameters,paths_dict,full_pathway,matplotlib_params,add_text_dict,add_text,windows))
        cleanup_button.bind('<ButtonRelease-1>',partial(tight_layout_on_release,windows,matplotlib_params))
        CreateToolTip(cleanup_button, text = 'Sometimes the y-label is cut off from the graph (does not apply to the final figure).\n'
                 'Tight Layout often fixes this issue, but the coordinate system displayed above is no longer consistent.\n'
                 'Previewing by holding this button will allow you to see how the labels will look like in the figure you will save,\n'
                 'without having to save it.')
        introduction = tk.Label(windows['graph_window'],text='''Hold down CTRL to check your cursor's coordinates.''')
        introduction.grid(row=1,column=1)
        destroy = update_graph(param_dict,parameters,paths_dict,full_pathway,windows['graph_window'],matplotlib_params,add_text_dict,add_text)
        windows['graph_window'].bind('<Control-Motion>',partial(mouse_coordinates,windows['graph_window'],execution,matplotlib_params['ax1'],param_dict))
        windows['graph_window'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'graph_window',matplotlib_params,'canvas'))
        if destroy:
            quit_window(windows,'graph_window',matplotlib_params,'canvas')
    return destroy

def add_text_converter(add_text_dict,add_text):
    dict_sorted = dict(sorted(add_text_dict.items()))
    add_text.clear()
    for index,value in dict_sorted.items():
       appends = [float(value[0]),float(value[1])]
       for params in value[2:]:
          appends.append(params)
       add_text.append(appends)


def generate_figure(param_dict,parameters,paths_dict,full_pathway,add_text_dict,add_text,save=False):
    make_parameters(param_dict,parameters)
    add_text_converter(add_text_dict,add_text)
    list_of_dicts,labelling = convert_path_to_list(paths_dict,full_pathway)
    pickle.dump([parameters,full_pathway,list_of_dicts,labelling,add_text],open( ".cache.fancy", "wb" ))
    if not list_of_dicts:
        messagebox.showerror('Error!',f"No pathway is defined, please define at least one.")
        ax1,fig = None,None
    else:
        if save:
          old_params = []
          for i in parameters:
              if i == 'visual=True':
                  old_params.append('visual=False')
              else:
                  old_params.append(i)
          parameters = old_params
          
        ax1,fig,savename = fancy_plots.init('Gibbs',parameters,full_pathway,list_of_dicts,labelling=labelling,add_text=add_text)
        if save:
            messagebox.showinfo('Save',f"Figure '{savename}' saved successfully.")
    
    return ax1,fig

def window(root):
    root['root'] = tk.Tk()
    root['root'].geometry('835x350')
    root['root'].title('Fancy Plots GUI')

def tkdestroy(window):
    try:
        window.destroy()
    except tk.TclError:
         pass

def quit_window(root,keyword,additional_dict=None,to_delete=None):
    tkdestroy(root[keyword])
    root[keyword] = None
    if additional_dict:
       del additional_dict[to_delete] 

def destroyall(root):
    for name,value in root.items():
      if value:
        tkdestroy(value)

def ok_click(root,name,location,event=None):
    filename=name.get()
    if not os.path.exists(f'{location}/{filename}.fancy'):
        os.system(f'mv {location}/.cache.fancy {location}/{filename}.fancy')
        destroyall(root)
    else:
        overwrite = messagebox.askquestion('Overwrite',f"Are you sure you want to overwrite {filename}.fancy?")
        if overwrite == 'yes':
          os.system(f'mv {location}/.cache.fancy {location}/{filename}.fancy')
          destroyall(root)
        else:
          pass

def quitall(root,location):
    save_cache = None
    if os.path.exists(f"{location}/.cache.fancy"):
        save_cache = messagebox.askyesnocancel('Quit',f"Do you want to save all settings to a file?")
        if save_cache and os.path.exists(f"{location}/.cache.fancy"):
            if 'savename'  in root:
                tkdestroy(root['savename'])
            root['savename'] = tk.Tk()
            root['savename'].geometry('300x75')
            root['savename'].title('Save cached information')

            label_name =tk.Label(root['savename'],text='Name of the file: ').grid(row=1,column=1)
            name       =tk.Entry(root['savename'],width=20,justify=tk.CENTER)
            name.grid(row=1,column=2)
            name_button= tk.Button(root['savename'],text='Save File',command=partial(ok_click,root,name,location)).grid(row=2,column=1,columnspan=2)

        elif save_cache is False:
           os.system(f'rm {location}/.cache.fancy')
           destroyall(root)

        else:
           pass
    
    else:
       destroyall(root)

def entry_tk(root,row,column,dictionary,name,message):
    dictionary[f'{name}_label']         = tk.Label(root,text=name)
    dictionary[f'{name}_label_row']     = row
    dictionary[f'{name}_label_column']  = column
    dictionary[name]                    = tk.Entry(root,width=10,justify=tk.CENTER)
    dictionary[f'{name}_row']           = row
    dictionary[f'{name}_column']        = column+1
    dictionary[f'{name}_label_message'] = message

def pack_entries(dictionary):
    for name, value in dictionary.items():
        if name.split('_')[-1] == 'label' or name.split('_')[-1] == 'loc' or name.split('_')[-1] == 'dec' or name.split('_')[-1] == 'min' or name.split('_')[-1] == 'double' or len(name.split('_')) == 1:
            value.grid(row=dictionary[f'{name}_row'],column=dictionary[f'{name}_column'])
        if name.split('_')[-1] == 'label':
            CreateToolTip(value, text = dictionary[f'{name}_message'])
   

def make_parameters(dictionary,parameters):
    parameters.clear()
    parameters.append('visual=True')
    for name, value in dictionary.items():
        if name.split('_')[-1] == 'loc' or name.split('_')[-1] == 'dec' or name.split('_')[-1] == 'min' or name.split('_')[-1] == 'double' or len(name.split('_')) == 1:
            if (name == 'xlabel' or name == 'ylabel' or name== 'title') and value.get():
                #parameters.append(f"{name}={value.get()}")
                parameters.append(f"{name}={value.get().encode().decode('unicode_escape')}")
            elif (name == 'xlabel' or name == 'ylabel') and not value.get():
                parameters.append(f"{name}={' '}")
            elif not value.get():
                pass
            else:        
                parameters.append(f'{name}={value.get()}')

def default_settings(dictionary,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text):
    for name, value in dictionary.items():
      if name.split('_')[-1] == 'loc' or name.split('_')[-1] == 'dec' or name.split('_')[-1] == 'min' or name.split('_')[-1] == 'double' or len(name.split('_')) == 1:
        dictionary[name].delete(0,tk.END)
    dictionary['font'].insert(0,'Sans Serif')
    dictionary['fontsize'].insert(0,'14')
    dictionary['linewidth'].insert(0,'2.25')
    dictionary['savename'].insert(0,'test.pdf')
    dictionary['dpi'].insert(0,'1200')
    dictionary['tick_loc'].insert(0,'out')
    dictionary['tick_dec'].insert(0,'2')
    dictionary['tick_min'].insert(0,'1')
    dictionary['tick_double'].insert(0,'False')
    dictionary['colors'].insert(0,'maroon,dodgerblue,m,g,k')
    dictionary['legend_loc'].insert(0,'best')
    dictionary['xlabel'].insert(0,'Reaction Coordinate')
    dictionary['ylabel'].insert(0,r'$\Delta$G (eV)')
    dictionary['boxsize'].insert(0,r'6,4.5')
    make_parameters(dictionary,parameters)
    if windows['graph_window']:
      update_graph(dictionary,parameters,paths_dict,full_pathway,windows['graph_window'],matplotlib_params,add_text_dict,add_text)
      execution['initial_coordinates'] = tk.Label(windows['graph_window'],text=f'x=0.000 y=0.000')
      execution['initial_coordinates'].grid(row=10,column=1,sticky=tk.S)
      windows['graph_window'].bind('<Control-Motion>',partial(mouse_coordinates,windows['graph_window'],execution,matplotlib_params['ax1'],dictionary))
      scale_window_x,scale_window_y=obtain_boxsizes(parameters)
      windows['graph_window'].geometry(f'{int((120*scale_window_x)+80)}x{int((120*scale_window_y)+80)}')

def pickle_settings(dictionary,parameters):
    for name, value in dictionary.items():
      if name.split('_')[-1] == 'loc' or name.split('_')[-1] == 'dec' or name.split('_')[-1] == 'min' or name.split('_')[-1] == 'double' or len(name.split('_')) == 1:
        dictionary[name].delete(0,tk.END)
    for name in parameters:
      name=name.split('=')
      if name[0] == 'visual':
        pass
      else: 
        dictionary[name[0]].insert(0,name[1])

def empty_line(root,row,column):
    tk.Label(root,text='').grid(row=row,column=column)

def get_path(path,path_names,event=None):
    path.clear()
    for i in path_names.get().split(','):
        string=''
        for j in i:
            if j==' ':
                pass
            else:
                string+=j
        path.append(string)

def save_values(pdict,def_entr,pathway,windows):
    success = True
    for name, value in def_entr.items():
        if name.split('_')[-1] == 'label' or name=='title':
            pass
        elif name=='legendlabel':
            #pdict[f'{pathway}-legendlabel'] = value.get()
            pdict[f'{pathway}-legendlabel'] = value.get().encode().decode('unicode_escape')
        else:
            try:
                pdict[f'{pathway}_{name}'] = float(value.get())
            except ValueError:
                if value.get():
                  success = False
                  messagebox.showwarning('Number not recognized!',f"'{value.get()}' is not recognized as a number, region '{name}' will be ignored for this pathway.")
                  windows['energy_window'].lift()
                  windows['energy_window'].after(1, lambda: windows['energy_window'].focus_force())
                if f'{pathway}_{name}' in pdict:
                   del pdict[f'{pathway}_{name}']
    if success:
        messagebox.showinfo('Save',f"Values on {pathway} saved successfully. Click on 'Show Graph' to create or update figure.")
        windows['energy_window'].lift()
        windows['energy_window'].after(1, lambda: windows['energy_window'].focus_force())

def make_path_windows(root,variable,full_pathway,paths_dict,define_entries,keyword=None):
    root['energy_window'].geometry(f'400x{160+(15*(len(full_pathway)//2))}')
    pathway = variable.get()

    empty_line(root['energy_window'],row=2,column=1)
    if len(full_pathway) < 3:
        define_entries['title'] = tk.Label(root['energy_window'],text=pathway).grid(row=3,column=1,columnspan=100,sticky=tk.W)
    else:
        define_entries['title'] = tk.Label(root['energy_window'],text=pathway)
        define_entries['title'].grid(row=3,column=1,columnspan=100)
        define_entries['title'].configure(anchor="center")
    xy = [4,-1]
    multiple=0
    for i in full_pathway:
        if multiple % 3 !=0 or multiple == 0:
            xy[1]+=2
        else:
            xy[0]+=1
        define_entries[f'{i}_label'] = tk.Label(root['energy_window'],text=i).grid(row=xy[0],column=xy[1])
        define_entries[i]            = tk.Entry(root['energy_window'],width=12)
        if paths_dict:
          try:
            define_entries[i].insert(0,str(float(paths_dict[f'{pathway}_{i}'])))
          except (ValueError, KeyError):
            pass
        define_entries[i].grid(row=xy[0],column=xy[1]+1)
        if xy[1]==5:
            xy[1]=1
        multiple+=1
    empty_line(root['energy_window'],row=xy[0]+1,column=1)
    tk.Label(root['energy_window'],text="Path Label \n (Legend)").grid(row=97,column=1,columnspan=3,sticky=tk.W)
    define_entries['legendlabel'] = tk.Entry(root['energy_window'],width=25)
    define_entries['legendlabel'].grid(row=97,column=1,columnspan=10,sticky=tk.W,padx=70)
    empty_line(root['energy_window'],row=98,column=1)
    empty_line(root['energy_window'],row=99,column=1)
    legendlabel = f'{pathway}-legendlabel'
    if legendlabel in paths_dict:
        #define_entries['legendlabel'].insert(0,str(paths_dict[f'{pathway}-legendlabel']))
        define_entries['legendlabel'].insert(0,str(paths_dict[legendlabel]).encode())
    save_button = tk.Button(root['energy_window'],text='Save',command=partial(save_values,paths_dict,define_entries,pathway,root)).grid(row=1,column=1,columnspan=10,sticky=tk.W,padx=105)

def define_free_energies(full_pathway,full_pathway_names,npathways,paths_dict,windows):
    if windows['energy_window']:
        quit_window(windows,'energy_window')
    if not full_pathway_names.get():
        messagebox.showerror("Divisions not found!","Please define the mechanism's divisions in 'Full Mechanism Divisions' entry box.")
        paths_dict.clear()
    else:
        windows['energy_window'] = tk.Toplevel()
        windows['energy_window'].geometry('350x80')
        windows['energy_window'].title('Fancy Plots - Energy Declaration')
        the_menu=make_textmenu(windows['energy_window'])
        windows['energy_window'].bind_class("Entry", "<Button-3><ButtonRelease-3>", partial(show_textmenu,the_menu))
        windows['energy_window'].bind_class("Entry", "<Control-q>", partial(callback_select_all,windows['root']))
        get_path(full_pathway,full_pathway_names)
        npaths = npathways.get()
        params = [f'Pathway {str(i+1)}' for i in range(int(npaths))]

        for string in list(paths_dict):
          number=string.split(" ")[1].split("_")[0].split("-")[0]
          if float(number) > int(npaths):
            del(paths_dict[string])

        define_entries = {}
        variable = tk.StringVar(windows['energy_window'])
        dropdown = tk.OptionMenu(windows['energy_window'], variable, *params,command=partial(make_path_windows,windows,variable,full_pathway,paths_dict,define_entries))
        variable.set('Pathway 1')
        make_path_windows(windows,variable,full_pathway,paths_dict,define_entries)

        dropdown.grid(row=1,column=1,columnspan=10,sticky=tk.W)

        windows['energy_window'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'energy_window'))

def make_parameters_and_plot(param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text,event=None):
    make_parameters(param_dict,parameters)
    if windows['graph_window']:
      update_graph(param_dict,parameters,paths_dict,full_pathway,windows['graph_window'],matplotlib_params,add_text_dict,add_text)
      execution['initial_coordinates'] = tk.Label(windows['graph_window'],text=f'x=0.000 y=0.000')
      execution['initial_coordinates'].grid(row=10,column=1,sticky=tk.S)
      windows['graph_window'].bind('<Control-Motion>',partial(mouse_coordinates,windows['graph_window'],execution,matplotlib_params['ax1'],param_dict))
      scale_window_x,scale_window_y=obtain_boxsizes(parameters)
      windows['graph_window'].geometry(f'{int((120*scale_window_x)+80)}x{int((120*scale_window_y)+80)}')

def open_image(name,xsize,ysize):
   location = os.path.realpath(__file__)
   split_location = location.split('/')
   del split_location[-1]
   split_location.extend(['fancy_plots_images',name])
   image_location = '/'.join(split_location)
   image = Image.open(image_location)
   image = image.resize((xsize,ysize),Image.ANTIALIAS)
   image = ImageTk.PhotoImage(image)
   return image

def mpl_palette(windows,images):
    if windows['matplotlib_palette']:
        quit_window(windows,'matplotlib_palette')

    windows['matplotlib_palette'] = tk.Toplevel()
    windows['matplotlib_palette'].geometry('553x800')
    windows['matplotlib_palette'].title('Fancy Plots - Color Palettes')

    '''color_scheme has to be added to a list so it's not garbage collected, and consequently not shown to the screen.
       Fixes to this issue are: Declare color_scheme as a global variable - not recommended. 
                                Create a class and declare color_scheme as self.color_scheme
                                Create a list and add this image to it. 
       Creating classes would simplify this code and get rid of most of lists and dictionaries, but that is something to be done in the future, 
       as the original fancy plots was written class-free and has to be rewritten almost entirely from scratch.'''

    color_scheme = open_image('color_palette.png',553,771)
    images[0]    = color_scheme
    mpl_image    = tk.Label(windows['matplotlib_palette'],image=images[0]).grid(row=0,column=0,columnspan=3)

    source    = tk.Label(windows['matplotlib_palette'],text='Source: matplotlib.org').grid(row=1,column=1)

    windows['matplotlib_palette'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'matplotlib_palette'))

def instructions(windows,images):
    if windows['instructions_window']:
        quit_window(windows,'instructions_window')

    windows['instructions_window'] = tk.Toplevel()
    windows['instructions_window'].geometry('1400x788')
    windows['instructions_window'].title('Fancy Plots - Instructions')
    tutorial          = open_image('fancy_plots_tutorial.png',1400,788)
    images[1]         = tutorial
    tutorial_image    = tk.Label(windows['instructions_window'],image=images[1]).grid(row=0,column=0,columnspan=3) 
    windows['instructions_window'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'instructions_window'))


def reorder_save(paths_dict,old,new,param_dict,parameters,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text,event=None):
    try:
        _new_arrangement = new.get().split(',')
        new_arrangement = []
        for i in _new_arrangement:
            int(i)
            if i not in old:
                raise ValueError
            elif i in new_arrangement:
                pass
            else:
                new_arrangement.append(i)
        if len(old) != len(new_arrangement):
            raise ValueError
           
    except ValueError:
        messagebox.showerror("Error!","Oh, no... Something went wrong. Make sure you type your desired order in this format --> 3,1,2 and that the path number is defined.")
        windows['reorder'].lift()
        windows['reorder'].after(1, lambda: windows['reorder'].focus_force())
        new_arrangement = False
  
    if new_arrangement:
        buffer = {}
        rearrange = {}
        for i,j in zip(old,new_arrangement):
            rearrange[j] = i
        for name in paths_dict:
            if name.split('_')[0].split(' ')[0] == 'Pathway' and not name.split('-')[-1] == 'legendlabel':
                number = name.split('_')[0].split(' ')[1]
                if name.split('_')[1].lower() == 'ts':
                   label = f"{name.split('_')[1]}_{name.split('_')[2]}"
                else:
                   label  = name.split('_')[1]
                buffer[f'Pathway {rearrange[number]}_{label}'] = paths_dict[f'Pathway {number}_{label}']
            elif name.split('-')[-1] == 'legendlabel':
                number = name.split('-')[0].split(' ')[1]
                buffer[f'Pathway {rearrange[number]}-legendlabel'] = paths_dict[f'Pathway {number}-legendlabel']
        paths_dict.clear()
        for name in buffer:
            paths_dict[name] = buffer[name] 
        graph(param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text)
        quit_window(windows,'reorder')        
def reorder_paths(paths_dict,windows,param_dict,parameters,full_pathway,matplotlib_params,execution,add_text_dict,add_text):
    if windows['reorder']:
        quit_window(windows,'reorder')
    npaths = 0
    namelist = []
    for name in paths_dict:
        if name.split('_')[0].split(' ')[0] == 'Pathway' and not name.split('-')[-1] == 'legendlabel':
            if int(name.split('_')[0].split(' ')[1]) >= npaths:
                npaths = int(name.split('_')[0].split(' ')[1])
            storename = f"Pathway {name.split('_')[0].split(' ')[1]}"
            if storename not in namelist:  
                namelist.append(storename)
    if npaths<2:
        messagebox.showerror("Is there a Pathway?","Please define at least 2 pathways.")
    else:
        windows['reorder'] = tk.Toplevel()
        windows['reorder'].geometry('250x100')
        windows['reorder'].title('Fancy Plots - Reordering pathways')
        windows['reorder'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'reorder'))
        old = []
        for i in namelist:
           old.append(i.split(' ')[1])
        old_arrangement_label=tk.Label(windows['reorder'],text='Old Arrangement').grid(row=1,column=1,sticky=tk.W)
        new_arrangement_label=tk.Label(windows['reorder'],text='New Arrangement').grid(row=1,column=2,sticky=tk.E,padx=30)
        old_arrangement=tk.Label(windows['reorder'],text=','.join(old)).grid(row=2,column=1,sticky=tk.W)
        new_arrangement=tk.Entry(windows['reorder'],width=10,justify=tk.CENTER)
        new_arrangement.grid(row=2,column=2)
        empty_line(windows['reorder'],row=3,column=1)
        press_reorder = tk.Button(windows['reorder'],text='Reorder',command=partial(reorder_save,paths_dict,old,new_arrangement,param_dict,parameters,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text)).grid(row=4,column=1,columnspan=2)
        windows['reorder'].bind('<Return>',partial(reorder_save,paths_dict,old,new_arrangement,param_dict,parameters,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text))       
        windows['reorder'].protocol('WM_DELETE_WINDOW', partial(quit_window,windows,'reorder'))

def show_text_def(number_text,add_text_dict,add_text_widget,add_text_x,add_text_y,add_text_color,add_text_fontsize,event=None):
   nth  = number_text.get()
   add_text_x.delete(0,tk.END)
   add_text_y.delete(0,tk.END)
   add_text_color.delete(0,tk.END)
   add_text_fontsize.delete(0,tk.END)
   add_text_widget.delete(0,tk.END)
   if add_text_dict:
      for index,value in add_text_dict.items():
         if int(index) == int(nth):
           add_text_x.insert(0,str(value[0]))
           add_text_y.insert(0,str(value[1]))
           #add_text_widget.insert(0,str(value[2]))
           add_text_widget.insert(0,str(value[2]).encode())
           if len(value) >= 4:
             add_text_color.insert(0,str(value[3]).split('=')[1])
           if len(value) == 5:
             add_text_fontsize.insert(0,str(value[4]).split('=')[1])

def add_text_def(number_text,add_text_dict,add_text_widget,add_text_x,add_text_y,add_text_color,add_text_fontsize,windows,param_dict,parameters,paths_dict,full_pathway,matplotlib_params):
   #text = add_text_widget.get()
   text = add_text_widget.get().encode().decode('unicode_escape') #recognizes \n
   x    = add_text_x.get()
   y    = add_text_y.get()
   nth  = number_text.get()
   if add_text_color.get():
      color= f'color={add_text_color.get()}'
   else:
      color='color=k'

   if add_text_fontsize.get():
      fontsize=f'fontsize={add_text_fontsize.get()}'
      add_text_dict[nth] = [x,y,text,color,fontsize]
   else:
      add_text_dict[nth] = [x,y,text,color]

   _continue          = True


   if not add_text_dict[nth][2]:
      del add_text_dict[nth]
   elif add_text_dict[nth][2] and (not isinstance(add_text_dict[nth][0],float) or not isinstance(add_text_dict[nth][1],float)):
      try:
         float(x)
         float(y)
         messagebox.showinfo("Save",f"Text {nth} saved successfully.")
      except ValueError:
         messagebox.showerror("Error!","Text is defined but at least one coordinate is not set or not valid. Please set your desired x AND y coordinates using numbers.")
         del add_text_dict[nth]
         _continue = False     
   if windows['graph_window'] and _continue:
      update_graph(param_dict,parameters,paths_dict,full_pathway,windows['graph_window'],matplotlib_params,add_text_dict,add_text)
      execution['initial_coordinates'] = tk.Label(windows['graph_window'],text=f'x=0.000 y=0.000')
      execution['initial_coordinates'].grid(row=10,column=1,sticky=tk.S)
      windows['graph_window'].bind('<Control-Motion>',partial(mouse_coordinates,windows['graph_window'],execution,matplotlib_params['ax1'],param_dict))
      scale_window_x,scale_window_y=obtain_boxsizes(parameters)
      windows['graph_window'].geometry(f'{int((120*scale_window_x)+80)}x{int((120*scale_window_y)+80)}')

def run():
    window(windows)
    here = os.getcwd()
    if len(sys.argv)>2:
      print('Error: More than 1 argument is not supported.')
      sys.exit(2)
    
    the_menu = make_textmenu(windows['root'])
    windows['root'].bind_class("Entry", "<Button-3><ButtonRelease-3>", partial(show_textmenu,the_menu)) #right mouse click popup cut, copy and paste menu
    windows['root'].bind_class("Entry", "<Control-q>", partial(callback_select_all,windows['root']))    # select all key binding - q
    windows['root'].protocol('WM_DELETE_WINDOW', partial(quitall,windows,here))          #redefine x button
    windows['energy_window']       = None
    windows['graph_window']        = None
    windows['reorder']             = None
    windows['tight_layout']        = None
    windows['instructions_window'] = None
    windows['matplotlib_palette']  = None
    matplotlib_params              = {}
    param_dict                     = {}
    paths_dict                     = {}
    add_text_dict                  = {}
    parameters                     = []
    add_text                       = []
    full_pathway                   = []
    images                         = [None,None]

    #Parameters
    param_title = tk.Label(windows['root'],text='Parameters',font='Helvetica 11 bold').grid(row=1,column=1, columnspan=30)
    empty_line(windows['root'],row=2,column=1)
    entry_tk(windows['root'],row=3,column=1,dictionary=param_dict,name='boxsize',message='Defines the width:height ratio of the graph by specifying *width*,*height* in inches.')  
    entry_tk(windows['root'],row=3,column=3,dictionary=param_dict,name='xaxis',message='Defines a non-default x-range by specifying *xmin*,*xmax*.')  
    entry_tk(windows['root'],row=3,column=5,dictionary=param_dict,name='yaxis',message='Defines a non-default y-range by specifying *ymin*,*ymax*.')    
    entry_tk(windows['root'],row=3,column=7,dictionary=param_dict,name='font',message="Defines a non-default font type. Currently supported: 'serif', 'sans-serif', 'monospace', 'cursive', 'fantasy'.")
    entry_tk(windows['root'],row=4,column=1,dictionary=param_dict,name='xlabel',message='Defines a label for the x-axis. If no label is desired, leave this space in blank.') 
    entry_tk(windows['root'],row=4,column=3,dictionary=param_dict,name='ylabel',message='Defines a label for the y-axis. If no label is desired, leave this space in blank.')
    entry_tk(windows['root'],row=4,column=5,dictionary=param_dict,name='tick_loc',message="Defines the tick's location - inside or outside of the graph's margin. Supports 'in' and 'out' keywords.")
    entry_tk(windows['root'],row=4,column=7,dictionary=param_dict,name='tick_dec',message="Defines the tick's decimal numbers. If '2' is stated, ticks will show e.g. 1.00, 2.00,...")
    entry_tk(windows['root'],row=5,column=1,dictionary=param_dict,name='tick_min',message='Defines how many minor ticks between the major ones are desired.')
    entry_tk(windows['root'],row=5,column=3,dictionary=param_dict,name='tick_double',message='If True, ticks will be shown on the right-hand side as well.')    
    entry_tk(windows['root'],row=5,column=5,dictionary=param_dict,name='fontsize',message="Size of the font for x-axis, y-axis and paths' labels. Title will have a value of the fontize+1 and additional text fontsize-2.") 
    entry_tk(windows['root'],row=5,column=7,dictionary=param_dict,name='linewidth',message="Width of the paths' lines, graph's margins are affected as well.") 
    entry_tk(windows['root'],row=6,column=1,dictionary=param_dict,name='legend_loc',message='Location of the legends. Accepted keywords: best, upper right, upper left, lower right, lower left, upper center, lower center, center left, center right, center.')     
    entry_tk(windows['root'],row=6,column=3,dictionary=param_dict,name='colors',message="Color palettes are shown if 'Color Palettes' button is clicked.")     
    entry_tk(windows['root'],row=6,column=5,dictionary=param_dict,name='title',message='States the title of the graph, leave this space in black if no title is desired.')   
    entry_tk(windows['root'],row=6,column=7,dictionary=param_dict,name='savename',message='Saves the figure with its corresponding extension (png,jpg,pdf,...).')       
    entry_tk(windows['root'],row=7,column=1,dictionary=param_dict,name='dpi',message='Dots per inch defines the resolution of the figure. This number will do nothing for pdf, svg and eps formats since these are vector images.')
    entry_tk(windows['root'],row=7,column=3,dictionary=param_dict,name='visual',message='Visual is enabled by default and cannot be changed. This keyword is recognized by fancy plots to enable GUI.')   
    pack_entries(param_dict)

    param_dict['visual'].insert(0,'True')
    param_dict['visual'].config(state='disabled')
    npaths = tk.IntVar(windows['root'])
    empty_line(windows['root'],row=8,column=1)
    full_pathway_names_label = tk.Label(windows['root'],text="Full Mechanism Divisions: ").grid(row=9,column=1,columnspan=3)
    npaths_label             = tk.Label(windows['root'],text="Number of Paths: ").grid(row=9,column=4,columnspan=1)
    full_pathway_names       = tk.Entry(windows['root'],width=30,justify=tk.CENTER)
    full_pathway_names.grid(row=10,column=1,columnspan=3)
    npathways = tk.Spinbox(windows['root'],from_=1,to=30,textvariable=npaths,width=3)
    npathways.grid(row=10,column=4,padx=0)
    npaths.set(1)

    energies                 = tk.Button(windows['root'], text='Define Gibbs Free Energies', command=partial(define_free_energies,full_pathway,full_pathway_names,npaths,paths_dict,windows)).grid(row=10,column=5,columnspan=2,padx=10)

    reorder                  = tk.Button(windows['root'], text='Reorder Pathways', command=partial(reorder_paths,paths_dict,windows,param_dict,parameters,full_pathway,matplotlib_params,execution,add_text_dict,add_text)).grid(row=10,column=5,columnspan=9,padx=196)
    default_settings(param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text)
    
    for name, value in param_dict.items():
        if name.split('_')[-1] == 'loc' or name.split('_')[-1] == 'dec' or name.split('_')[-1] == 'min' or name.split('_')[-1] == 'double' or len(name.split('_')) == 1:
            value.bind('<Return>',partial(make_parameters_and_plot,param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text))
    empty_line(windows['root'],row=11,column=1)

    number_text           = tk.IntVar(windows['root'])
    add_text_label        = tk.Label(windows['root'],text="Additional Text: ").grid(row=12,column=1,columnspan=3)
    add_text_x_label      = tk.Label(windows['root'],text="X Coordinate: ").grid(row=12,column=4,columnspan=1)
    add_text_y_label      = tk.Label(windows['root'],text="Y Coordinate: ").grid(row=12,column=5,columnspan=1)
    add_text_color_label  = tk.Label(windows['root'],text="Color: ").grid(row=12,column=6,columnspan=1) 
    add_text_fontsize     = tk.Label(windows['root'],text="Fontsize: ").grid(row=12,column=7,columnspan=1)
    add_text_widget       = tk.Entry(windows['root'],width=30,justify=tk.CENTER)
    add_text_x            = tk.Entry(windows['root'],width=10,justify=tk.CENTER)
    add_text_y            = tk.Entry(windows['root'],width=10,justify=tk.CENTER)
    add_text_color        = tk.Entry(windows['root'],width=10,justify=tk.CENTER)
    add_text_fontsize     = tk.Entry(windows['root'],width=3,justify=tk.CENTER)
    n_text_label          = tk.Label(windows['root'],text="Text Number: ").grid(row=12,column=8,columnspan=1)
    n_text                = tk.Spinbox(windows['root'],from_=1,to=100,textvariable=number_text,width=3,command=partial(show_text_def,number_text,add_text_dict,add_text_widget,add_text_x,add_text_y,add_text_color,add_text_fontsize))
    n_text.bind('<Return>',partial(show_text_def,number_text,add_text_dict,add_text_widget,add_text_x,add_text_y,add_text_color,add_text_fontsize))
    save_text        = tk.Button(windows['root'],text='Save Text', command = partial(add_text_def,number_text,add_text_dict,add_text_widget,add_text_x,add_text_y,add_text_color,add_text_fontsize,windows,param_dict,parameters,paths_dict,full_pathway,matplotlib_params)).grid(row=13,column=9,columnspan=1,sticky=tk.W)
    add_text_widget.grid(row=13,column=1,columnspan=3)
    add_text_x.grid(row=13,column=4)
    add_text_y.grid(row=13,column=5)
    add_text_color.grid(row=13,column=6)
    add_text_fontsize.grid(row=13,column=7)
    n_text.grid(row=13,column=8)
    
    number_text.set(1)

    defaults     = tk.Button(windows['root'], text='Default Settings', command=partial(default_settings,param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text)).grid(row=14,column=1,columnspan=2,pady=30)
    create_graph = tk.Button(windows['root'], text='Show Graph', command=partial(graph,param_dict,parameters,paths_dict,full_pathway,windows,matplotlib_params,execution,add_text_dict,add_text)).grid(row=14,column=2,columnspan=10,pady=30,padx=70,sticky=tk.W)
    save_graph   = tk.Button(windows['root'], text='Save Graph', command=partial(generate_figure,param_dict,parameters,paths_dict,full_pathway,add_text_dict,add_text,save=True)).grid(row=14,column=3,columnspan=10,pady=30,padx=77,sticky=tk.W)
   

    matplotlib_palette_button = tk.Button(windows['root'],text='Color Palettes',command=partial(mpl_palette,windows,images)).grid(row=14,column=4,columnspan=10,pady=30,padx=94,sticky=tk.W)
    instructions_button = tk.Button(windows['root'],text='Instructions',command=partial(instructions,windows,images)).grid(row=14,column=5,columnspan=10,pady=30,padx=89,sticky=tk.W)
    
    if len(sys.argv)>1:
        cache=sys.argv[1]
        if os.path.exists(f'{here}/{cache}'):
            execute = True
            try:
              all_unpacked  = pickle.load(open(cache, "rb" ))
            except:
              messagebox.showwarning('Cache Warning!',f"File '{cache}' could not be opened - It either is not a .fancy file or it was damaged. It will be ignored.")
              execute=False
            if execute:
              if len(all_unpacked) != 5:
                messagebox.showwarning('Cache Warning!',f"File '{cache}' has wrong format. It will be ignored.")
              elif execute:
                pickle_parameters    = all_unpacked[0]
                pickle_full_pathway  = all_unpacked[1]
                pickle_list_dicts = all_unpacked[2]
                pickle_labelling     = all_unpacked[3]
                pickle_text      = all_unpacked[4]
                pickle_settings(param_dict,pickle_parameters)
                for index,value in enumerate(pickle_text):
                  text_params = [text_param for text_param in value]
                  add_text_dict[index+1] = text_params
                show_text_def(number_text,add_text_dict,add_text_widget,add_text_x,add_text_y,add_text_color,add_text_fontsize)
                npaths.set(len(pickle_list_dicts))
                full_pathway_names.insert(0,','.join(pickle_full_pathway)) #.encode shows characters
                for index, dictionary in enumerate(pickle_list_dicts):
                  for name,value in dictionary.items():
                    paths_dict [f'Pathway {index+1}_{name}'] = value
                  paths_dict [f'Pathway {index+1}-legendlabel'] = pickle_labelling[index]
                define_free_energies(full_pathway,full_pathway_names,npaths,paths_dict,windows)

    windows['root'].mainloop()        

if __name__ == '__main__' :
    run()

