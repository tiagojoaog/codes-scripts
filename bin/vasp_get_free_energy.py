#!/usr/bin/env python

import sys
from ase.io import read
from ase.thermochemistry import HarmonicThermo,IdealGasThermo
import os

vibrational_analysis=open('vib.txt', 'r')
v=vibrational_analysis.readlines()
i=1
freq=[]
GAS="false"
linear="false"
temp=298.15
pbar=1
spin_entry=0
symm_num=1
thresh=12
print(" ")
for trial in sys.argv:
        if trial=="gas":
                GAS="true"
                print("  Gas-phase approximation: activated  ")
                print("")
                print(" ---------------------------------------------------------------------------------------------------------- ")
                print(" GAS PHASE: 3 trans 3 rot (2rot if linear molecule) discarded, used ideal gas approx for contributions.")
                print("")
                print(" ---------------------------------------------------------------------------------------------------------- ")
                print("")
        elif trial=="linear":
                i+=1
                linear="true"
                print("  Linear-geometry: activated  ")
                print("")
        elif trial=="TS":
                i-=1
                print("  Transition-state: activated  ")
                print("")
        elif trial=="solv":
                thresh=100
                print("  Solution: activated  ")
                print("")
                print(" ---------------------------------------------------------------------------------------------------------- ")
                print(" In Solution: low vibrations are shifted to 100cm-1 instead.")
                print("")
                print(" ---------------------------------------------------------------------------------------------------------- ")
                print("")
        elif trial.split("=")[0]=="pressure":
                pbar=float(trial.split("=")[1])
        elif trial.split("=")[0]=="spin":
                print("You activated non-null spin.")
                print("Note: Only relevant for the gas-phase approximation.")
                print("")
                spin_entry=float(trial.split("=")[1])
                print("Spin is ",spin_entry)
                print("")
        elif trial.split("=")[0]=="temperature":
                temp=float(trial.split("=")[1])
        elif trial.split("=")[0]=="symmetry":
                print("You have stated a symmetry number, please make sure it corresponds to the molecule in question.")
                print("Note: Only relevant for the gas-phase approximation.")
                symm_num=float(trial.split("=")[1])
                print("")
                print("Symmetry number stated is ",symm_num)
                print("")
if "gas" not in sys.argv:
        print("")
        print(" ---------------------------------------------------------------------------------------------------------- ")
        print(" Harmonic limit: Harmonic energies of the adsorbate (all frequencies below "+str(thresh)+" cm^-1 are treated as "+str(thresh)+" cm^-1).")
        print("")
        print(" ---------------------------------------------------------------------------------------------------------- ")
        print(" ")
print(" ---------------------------------------------------------------------------------------------------------- ")
freq_cons=open("gibbs_freq_used.txt","w")
freq_cons.write("------- Frequencies used for Gibbs Free energy calculation ---------"+"\n")
freq_cons.write("\n")
for line in v:
        try:
                int(line.split()[0])
                i+=1
                if i>7 and GAS=="true":
                        word=line.split()
                        ""
                        " GAS PHASE: 3 trans 3 rot (2rot if linear molecule) discarded, used ideal gas approx for contributions."
                        ""
                        try:
                                var=float(word[2])
                                freq_cons.write(str(var)+"\n")
                                freq.append(var*0.0001239843)
                        except ValueError:
                                print("")
                                print(" Warning: Imaginary frequencies besides trans, rot or TS. Those were not taken into account, try to reoptimize the structure.")
                                print("")
                elif i>1 and GAS=="false":
                        word=line.split()
                        ""
                        " Harmonic limit: Harmonic energies of the adsorbate (all frequencies below "+str(thresh)+" cm^-1 are treated as "+str(thresh)+" cm^-1)."
                        ""
                        try:
                                var=float(word[2])
                                if float(word[2])< thresh:
                                        var=thresh
				#freq_cons.write(str(var)+"\n")
				#freq.append(var*0.0001239843)
                        except ValueError:
                                print("")
                                print(" Warning: Imaginary frequencies besides TS. Changed to "+str(thresh)+" cm^-1.")
                                print("")
                                var=thresh
                        freq_cons.write(str(var)+"\n")
                        freq.append(var*0.0001239843)
        except ValueError:
                continue
print(" ---------------------------------------------------------------------------------------------------------- ")
print("")
freq_cons.write("\n")
freq_cons.write("--------------------------------------------------------------------"+"\n")
freq_cons.write("--------------------------------------------------------------------"+"\n")
freq_cons.close()
vibrational_analysis.close()
#pot_energ=open('final.e', 'r')
#poten=float(pot_energ.read())
#print poten
#pot_energ.close()
atoms=read('in.traj')
#print freq

if GAS=="true" and linear=="false":
        thermo = IdealGasThermo(atoms=atoms,
                                vib_energies=freq,
                                geometry='nonlinear',
                                potentialenergy=0.0,
                                symmetrynumber=symm_num,spin=spin_entry)

        dg_ase = thermo.get_gibbs_energy(temp, pbar*100000, verbose=True)

elif GAS=="true" and linear=="true":
        thermo = IdealGasThermo(atoms=atoms,
                                vib_energies=freq,
                                geometry='linear',
                                potentialenergy=0.0,
                                symmetrynumber=symm_num,spin=spin_entry)

        dg_ase = thermo.get_gibbs_energy(temp, pbar*100000, verbose=True)

else:
        thermo = HarmonicThermo(vib_energies=freq, potentialenergy=0.0)

        dg_ase  = thermo.get_helmholtz_energy(temp,verbose=True)

f=open("dg_ase.log","w")
f.write(str(dg_ase))
f.close()
