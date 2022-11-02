#!/usr/bin/env python3

import os, sys
import numpy as np
import matplotlib.pyplot as plt

def read_d_occ(filename):
    i=-1
    with open(filename) as f:
        for line in f:
            if 'total' in line.split(' ') and 'charge' in line.split(' '):
                i+=1   
            elif i>-1:
                i+=1

            if i==4:
                return float(line.split(' ')[-4])



d_occ_gs = read_d_occ('calc0/OUTCAR')       #groundstate, folder names may be different

plot_U = np.arange(-0.2,0.2,0.05) # change this in case different intervals have been used

plot_d_bare  = np.empty(plot_U.size, dtype=float)
plot_d_inter = np.empty(plot_U.size, dtype=float)
plot_d_gs    = np.ones(plot_U.size, dtype=float)*d_occ_gs

for index,value in enumerate(plot_U):
    if abs(value) > 1e-4:
        bare  = 'calc'+str(round(value,2))+'bare'  #folder names may be different
        inter = 'calc'+str(round(value,2))+'inter'

        plot_d_bare[index]  = read_d_occ(f'{bare}/OUTCAR')    #non SCF
        plot_d_inter[index] = read_d_occ(f'{inter}/OUTCAR')  #SCF
    else:
        plot_d_bare[index]  = d_occ_gs
        plot_d_inter[index] = d_occ_gs


slope_b     = np.polyfit(plot_U,plot_d_bare,1)[0]
intercept_b = np.polyfit(plot_U,plot_d_bare,1)[1]
slope_i     = np.polyfit(plot_U,plot_d_inter,1)[0]
intercept_i = np.polyfit(plot_U,plot_d_inter,1)[1]

print('\n',plot_U)
print(plot_d_bare)
print(plot_d_inter,'\n')
print(f'Non SCF f(x)={slope_b:.2f}*x + {intercept_b:.2f}')
print(f'SCF     f(x)={slope_i:.2f}*x + {intercept_i:.2f}')

plt.plot(plot_U, plot_d_bare, ls='--', marker='o', c='r', label='non SCF')
plt.plot(plot_U, plot_d_inter, ls='--', marker='o', c='b', label='SCF')
plt.legend(loc='best')
plt.xlabel('U (eV)')
plt.ylabel('d occupancy')
plt.text(0.00, 6.00, f'non SCF f(x)={slope_b:.2f}*x + {intercept_b:.2f}')
plt.text(0.00, 5.95, f'SCF       f(x)={slope_i:.2f}*x + {intercept_i:.2f}')
plt.show()

print('\n')
print(f'U value is :{((1/slope_i)-(1/slope_b)):.2f} eV')
