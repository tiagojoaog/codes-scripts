#!/usr/bin/env python

#from matplotlib import pyplot as plt
import sys
import ase
#from ase.io import read
import numpy as np
#from ase.thermochemistry import HarmonicThermo,IdealGasThermo
#from turbomole_python import tm_get_gibbs_from_fil
import os
from ase.io import *

os.system("chgsum.pl AECCAR0 AECCAR2")
os.system("bader CHGCAR -ref CHGCAR_sum")
os.system("bader_analysis.py")
os.system("cat CHARGE_ANALYSIS.txt | awk '{print $7}'>charge.tmp")

while True:
  a=input("Check charge of molecule (input indices separated by comma e.g. 1,2,3): ")
  summation=0
  print(a.split(","))
  for i in a.split(","):
    j=0
    f=open("charge.tmp","r")
    for l in f:
      #print int(i)+1, j
      if j==int(i)+1:
        #print "success"
        summation=summation+float(l)
      j+=1
    f.close()
  print("")
  print("")
  print("Charge of molecule with indices "+a+" is "+str(summation))
  print("")
  print("")
  again=input("Again?(yes/no) ")
  print("")
  print("")
  if again=="yes":
    pass
  else:
    break

