#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os,sys


def Ulim(Urhe,OH,OOH,O2):
    for i in Urhe:
        O=2*OH+0.28
        ind_OH=OH-1*i
        ind_O=O-2*i-ind_OH
        ind_OOH=OOH-3*i-(O-2*i)
        ind_O2= O2-4*i-(OOH-3*i)
        global_O2 = O2-4*i

        dummy_list = [ind_OH, ind_O, ind_OOH, ind_O2, global_O2]
        Ulim_check = max(dummy_list)
        if abs(Ulim_check) < 1E-2:
            Ulim=Ulim_check
            return i

Ulim_list = []

Urhe = np.linspace(0,5,300)
OH_1d = np.linspace(0,3,300)
OOH_1d = np.linspace(3,5,300)
O2 = 4.92

OOH=[]
OH=[]

for indy,numy in enumerate(OOH_1d):
    dumx=[]
    dumy=[]
    ulim_dum=[]
    for indx,numx in enumerate(OH_1d):
        dumx.append(numx)
        dumy.append(numy)
        ulim_dum.append(Ulim(Urhe,numx,numy,O2))
    OH.append(dumx)
    OOH.append(dumy)
    Ulim_list.append(ulim_dum)


OH_line =np.linspace(0.2,1.6,10)
line_corr = [i+3.2 for i in OH_line]

fig, ax = plt.subplots(1, 1)

im=ax.pcolormesh(OH, OOH, Ulim_list, cmap='RdYlBu',vmin=1.2, vmax=3.6, shading='auto')
cbar=fig.colorbar(im, ax=ax)

cbar.set_label('Ulim (V)')

plt.plot(OH_line,line_corr,color='k', label='E(OOH)=E(OH*)+3.2 eV')

plt.xlabel('OH* (eV)')
plt.ylabel('OOH* (eV)')
plt.legend(loc='best')
plt.xlim([0,3])
plt.ylim([3,5])
plt.title('Heatmap for OER')

plt.show()

