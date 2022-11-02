#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.lines
import sys


def add_text_instances(add_text,visual):
    fontsize_label=fontsize-2
    if not add_text:
      pass
    else:
      for text in add_text:
        sett_count=0
        text_color='k'
        for setting in text:
          sett_count+=1
          if sett_count<=3:
            pass
          elif 'color' in setting.split('=')[0] or 'c' in setting.split('=')[0]:
            text_color=setting.split('=')[1]
          elif 'fontsize' in setting.split('=')[0]:
            fontsize_label=setting.split('=')[1]
        if not visual:
          plt.text(text[0],text[1],text[2],color=text_color,fontsize=fontsize_label,fontfamily=font)
        else:
          ax1.text(text[0],text[1],text[2],color=text_color,fontsize=fontsize_label,fontfamily=font)
        #plt.text(text[0],text[1],text[2],color=text_color,fontsize=fontsize,fontfamily=font)

def interpol(x0,y0):
  xx=list(np.linspace(x0[0],x0[1],50))
  prefac=y0[1]-y0[0]
  yy=[y0[0] + prefac *np.sin(0.5*np.pi*(x-x0[0])/(x0[1]-x0[0])) for x in xx ]
  xx2= list(np.linspace(x0[1],x0[2],50))
  prefac=y0[2]-y0[1]
  yy+=[y0[2] - prefac *np.sin(0.5*np.pi*(x-x0[0])/(x0[2]-x0[1])) for x in xx2 ]
  xx+=xx2
  return xx,yy

def create_list(mol,lists,color,zorder,label):
    count=-1
    mol_list_num=-1
    x_dash=[]
    y_dash=[]
    for i in mol:
        if i not in lists:
          print('\nWARNING: ',i,' not in at least 1 of the dictionaries.\n')
          mol_list_num+=1
          if 'TS'==i.split('_')[0].upper():
            count+=1
          else:
            count+=2
        elif 'TS'==i.split('_')[0].upper():
            x=[]
            y=[]
            count+=1
            #x.append(count)
            #y.append(lists[i])
            mol_list_num+=1
            index_check_before=1
            index_check_after=1
            index_before=mol_list_num-index_check_before
            index_after=mol_list_num+index_check_after
            trigger=0
            trigb=0
            triga=0
            while trigger==0:
              if mol[index_before] not in lists: 
                index_check_before+=1
                index_before=mol_list_num-index_check_before
                if index_before<0:
                  print("\nError! Please check if your pathway makes sense.")
                  print("       Structure before your transition state not found")
                  print("       for at least one of your pathways.\n")
                  sys.exit(2)
              else:
                trigb=1
              if mol[index_after] not in lists:
                index_check_after+=1
                index_after=mol_list_num+index_check_after
                if index_after>len(mol)-1:
                  print()
                  print("\nError! Please check if your pathway makes sense.")
                  print("       Structure after your transition state not found")
                  print("       for at least one of your pathways.\n")
                  sys.exit(3)
              else:
                triga=1
              if trigb==1 and triga==1:
                trigger=1
            a_cut=[count-0.99,count,count+0.99]
            b_cut=[lists[mol[index_before]],lists[mol[mol_list_num]],lists[mol[index_after]]]
            xx,yy = interpol(a_cut, b_cut)
            x.extend(xx)
            y.extend(yy)
            x_dash.extend(xx)
            y_dash.extend(yy)
            plot_cont(x,y,color,zorder,label)
        else:
            x=[]
            y=[]
            count+=1
            x.append(count)
            x_dash.append(count)
            count+=1
            x.append(count)
            x_dash.append(count)
            y.extend([lists[i],lists[i]])
            y_dash.extend([lists[i],lists[i]])
            plot_cont(x,y,color,zorder,label)
            if label:
              label=False
            else:    
              pass
            mol_list_num+=1
    plot_dash(x_dash,y_dash,color,zorder)
    return count

def plot_dash(x,y,color,zorder):
  dashedline=linewidth*0.7
  ax1.plot(x,y,'--',color=color,linewidth=dashedline,zorder=zorder)
def plot_cont(x,y,color,zorder,label):
  if label:
    ax1.plot(x,y,'-',color=color,linewidth=linewidth,zorder=zorder,label=label)
  else:
    ax1.plot(x,y,'-',color=color,linewidth=linewidth,zorder=zorder)

def gibbs_plt(dict_list,labelling,mol_list,colors):
  iteration=0
  for structs in dict_list:
    try:
      label=labelling[iteration]
    except:
      label=''
    try:
      c=colors[iteration]
    except:
      c='k'
      print('WARNING: Number of colors defined are not sufficient for the pathways defined. Taking the default color (black).')
    count=create_list(mol_list,structs,c,10-iteration,label=label) 
    iteration+=1
  return count

def general(dict_list,labelling,colors):
  iteration=0
  for structs in dict_list:
    try:
      label=labelling[iteration]
    except:
      label=''
    if len(structs)<3:
      print('Error: Add the linestyle as the last element of each list. E.g. \'-\'')
      sys.exit(4)
    try:
      ax1.plot(structs[0],structs[1],structs[2],color=structs[3],linewidth=linewidth,markersize=markersize,zorder=10-iteration,label=label)
    except:
      ax1.plot(structs[0],structs[1],structs[2],color=colors[iteration],linewidth=linewidth,markersize=markersize,zorder=10-iteration,label=label)
    iteration+=1


def init(kind,parameters,mol_list,dict_list,labelling,add_text):
  global ax1,fig,linewidth,font,fontsize,markersize

  variables = {'legend_loc':'best' , 'font':'Arial', 'fontsize':9 ,
               'markersize':7 , 'linewidth': 1.5, 'markeredgewidth':0.3 ,
               'xlabel':'Reaction Coordinate', 'ylabel':u'$\Delta$G (eV)' ,
               'xscale':'', 'yscale':'', 'savename':'test.png' , 'dpi':1200,
               'title':'' , 'visual':False, 'tick_loc':'out', 'tick_dec':'2', 
               'tick_min':'1', 'tick_double':False, 'colors':["maroon","dodgerblue",'m','g','k']}

  changed_variables = []

  keys_string = ['font','title','legend_loc','savename', 'tick_loc', 'tick_dec', 'tick_min',
                     'xlabel','ylabel','xscale','yscale']

  keys_float  = ['fontsize','markersize','linewidth',
                                 'markeredgewidth','dpi']

  keys_ranges = ['xaxis','yaxis','boxsize']

  keys_booleans = ['visual','tick_double']

  if not parameters:
    print("\nWARNING: No parameters set. Using the default ones.\n")

  for name in parameters:

    keyword=name.split('=')
    var = keyword[0]
    string = keyword[1]

    if var in keys_string:
      globals()[var] = string
      changed_variables.append(var)
      print(var,':',string)
    elif var in keys_float:
      globals()[var] = float(string)
      changed_variables.append(var)
      print(var,':',string)
    elif var in keys_ranges:
      highly_specific_a=float(string.split(',')[0])
      highly_specific_b=float(string.split(',')[1])
      globals()[var+'1'] = highly_specific_a
      globals()[var+'2'] = highly_specific_b
      changed_variables.append(var)
      print(var+'1',':',highly_specific_a)
      print(var+'2',':',highly_specific_b)

    elif var in keys_booleans:
      if string == 'False':
        globals()[var] = False
      else:
        globals()[var] = True
      changed_variables.append(var)
      print(var,':',string)

    elif var == 'colors':
        colors = string.split(',')
        changed_variables.append(var)
        print(var,':',colors)
  
    else:
      
      try:
        print(var,': not set   ERROR: Keyword not recongized.')
      except:
        pass

  for checker in variables:
      if checker not in changed_variables:
          if checker in keys_string:
            globals()[checker] = variables[checker]
            print(checker,':',variables[checker])
          elif checker in keys_float:
            globals()[checker] = variables[checker]
            print(checker,':',variables[checker])
          elif checker in keys_booleans:
            globals()[checker] = variables[checker]
            print(checker,':',variables[checker])
          elif checker == 'colors':
            colors = variables['colors']
            print(checker,':',variables['colors'])
  print('\n \n')
  #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
  #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
  #plt.rcParams["font.family"] = font
  #plt.rcParams.update({'font.size': fontsize})
  #plt.rcParams.update({'legend.fontsize':fontsize})
  #plt.rcParams.update({'xtick.labelsize':fontsize})
  if not visual:
    print('\n \n \n')
    print("#######################################################################################################")
    print("#                                                                                                     #")
    print("#                                                                                                     #")
    print("#    *****    *****   ***    **   ****  **     **     ******   **       ***    **********  *******    #")
    print("#    **      **  **   **  *  **  *****   **  **       **   **  **      ** **   **********  **         #")
    print("#    ****   *** ***   **   ****  **        **         ** **    **     **   **      **      *******    #")
    print("#    **    **    **   **    ***  *****     **         **       ******  ** **       **           **    #")
    print("#    **   **     **   **     **   ****     **         **       ******   ***        **      *******    #")
    print("#                                                                                                     #")
    print("#    Version 0.1                                           Contact info:                              #")
    print("#                                                            -tiagojoaog@gmail.com                    #")
    print("#######################################################################################################")
    print('\n \n \n')
    plt.rcParams["font.family"] = font
    plt.rcParams.update({'font.size': fontsize})
    plt.rcParams.update({'legend.fontsize':fontsize})
    plt.rcParams.update({'xtick.labelsize':fontsize})
    fig, ax1 = plt.subplots()
    plt_plot = True
  else:
    from matplotlib import rc
    from matplotlib.figure import Figure
    matplotlib.use("TkAgg")

    if 'boxsize' in changed_variables:
        boxsize_tuple=(boxsize1,boxsize2)
    else:
        boxsize_tuple=(4,3)
                      
    fig     = Figure(figsize=boxsize_tuple,dpi=120)   #ideal 6,4

    rc('font',**{'family':font,'size':fontsize})
    rc('legend',**{'fontsize':fontsize})
    rc('xtick',**{'labelsize':fontsize})
    
    plt_plot = False

    ax1     = fig.add_subplot(111)

  if not visual:
    if 'boxsize' in changed_variables:
      fig.set_size_inches(boxsize1, boxsize2)
    else:
      fig.set_size_inches(4.0, 3.2)
  box_width=linewidth
  #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
  #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
  if kind.lower()=='gibbs':
    count=gibbs_plt(dict_list,labelling,mol_list,colors) 
  else:
    general(dict_list,labelling,colors)



  if 'xaxis' in changed_variables:
    ax1.set_xlim([xaxis1,xaxis2])
  else:
    if kind.lower()=='gibbs':
      ax1.set_xlim([-0.05,count+0.05])
  if 'yaxis' in changed_variables:
    ax1.set_ylim([yaxis1,yaxis2])
  for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(box_width)
    ax1.spines[axis].set_zorder(100)
  
  ax1.set_xlabel(xlabel, fontsize=fontsize,fontfamily=font)
  ax1.set_ylabel(ylabel, fontsize=fontsize,fontfamily=font)
  
  #ax1.xaxis.set_minor_locator(AutoMinorLocator(5))  #AutoMinorLocator(3) specifies 2 minor ticks between major ones
  ax1.yaxis.set_minor_locator(AutoMinorLocator(int(tick_min)+1))
  ax1.yaxis.set_major_formatter(FormatStrFormatter(f'%.{tick_dec}f'))
  ax1.tick_params(which='both', direction=tick_loc, labelsize=fontsize, width=linewidth)
  ax1.tick_params(which='major', length=4*linewidth)
  ax1.tick_params(which='minor', length=2.5*linewidth)
  if kind.lower()=='gibbs':
    ax1.xaxis.set_ticklabels([])
    ax1.xaxis.set_ticks([])
    #ax1.yaxis.set_ticklabels([])  #
    #ax1.yaxis.set_ticks([])       #
  else:
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))    
  #plt.xticks([373,423,473,523,573,623,673],fontsize=fontsize)
  #plt.yticks([0,20,40,60,80,100],fontsize=fontsize)
  if title:
    if plt_plot:
        plt.title(title,fontsize=fontsize+1,fontfamily=font)
    else:
        ax1.set_title(title,fontsize=fontsize+1,fontfamily=font)
  if xscale:
    ax1.set_xscale(xscale)
  if yscale:
    ax1.set_yscale(yscale)

  if tick_double:
    #ax1.twinx()
    ax2=ax1.secondary_yaxis('right')
    ax2.yaxis.set_minor_locator(AutoMinorLocator(int(tick_min)+1))
    ax2.yaxis.set_major_formatter(FormatStrFormatter(f'%.{tick_dec}f'))
    ax2.tick_params(which='both', direction='in', labelsize=fontsize, width=linewidth)
    ax2.tick_params(which='major', length=4*linewidth)
    ax2.tick_params(which='minor', length=2.5*linewidth)
    ax2.yaxis.set_ticklabels([])

  if plt_plot:
    plt.legend(loc=legend_loc,frameon=False,fontsize=fontsize-2)
  else:
    ax1.legend(loc=legend_loc,frameon=False,fontsize=fontsize-2)
  add_text_instances(add_text,visual)
  if not visual:
    if savename.split('.')[1].lower()=='png':
        plt.savefig(savename, bbox_inches="tight",dpi=dpi,transparent=True)
    else:
        plt.savefig(savename, bbox_inches="tight",transparent=True)
    plt.close()

  return ax1,fig,savename
    
