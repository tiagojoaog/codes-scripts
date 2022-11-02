#!/usr/bin/env python3

from classical_mechanics import *
from matplotlib import pyplot as plt

Analysis=orbit(frames=800,capture=True,draw=True)   #Use a value for frames to capture images and make a GIF, also enable capture
Analysis.analysis(traj_show_save='show',animate=True) #traj_show_save either shows or saves trajectory pdf and draw enables object animation
                                                      #Animates trajectories
#Analysis.clear_cache()                               #Deletes pickle (cache) files

Analysis.close()  #closes the program
