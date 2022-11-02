#!/usr/bin/env python3

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List
from matplotlib import pyplot as plt
import numpy as np
import os, sys

'''

This script can be generalized if
folder names are not hardcoded.

e.g. script is imported and takes
foldernames as a parameter in 
get_properties().run().

'''

class FUNCTIONAL(Enum):
    PBE    = auto()
    PBED3  = auto()

class SOFTWARE(Enum):
    QE     = auto()
    VASP   = auto()

class PP(Enum):
    LANNA  = auto()
    SAMIRA = auto()


@dataclass(order=True)
class structure():
    energy      : float
    name        : str
    folder_name : str
    molecule    : str
    Nb          : str
    software    : SOFTWARE
    ecut        : float
    functional  : FUNCTIONAL
    pseudopot   : PP

    def __post_init__(self) ->None:
        object.__setattr__(self, 'sort_index', self.energy)

    def __str__(self) -> str:
        return f"""
{self.name}:
  * Folder                 --> {self.folder_name}
  * Name                   --> {self.name}
  * Functional             --> {self.functional}
  * Software               --> {self.software}
  * Energy Cutoff (eV)     --> {self.ecut}
  * Pseudopotential        --> {self.pseudopot}
  * PP Nb                  --> {self.Nb}
  * Structure type         --> {self.molecule}
  * Electronic Energy (eV) --> {self.energy}
                """

class get_properties():
    def __init__(self):
        self._list     = []
        #self.QE_mols   = '/home/tiago.ferreiragoncal/Lanna_test/QE/molecules'
        #self.VASP      = '/home/tiago.ferreiragoncal/Lanna_test/VASP'
        #self.QE        = '/home/tiago.ferreiragoncal/Lanna_test/QE'

    def get_e(self,name):
        energy=None
        if os.path.exists(name):
            for e in open(name):
                try:
                    energy=float(e)
                except ValueError:
                    pass
        return energy

    def call_obj(self,*kwargs):
        for obj in self._list:
            if obj.name in kwargs and obj.functional in kwargs and obj.software in kwargs and obj.pseudopot in kwargs and obj.Nb in kwargs and obj.ecut in kwargs:
                return obj

    def iterate_folders(self,add_name='',*,directory,func,soft,pp,molecule,Nb,ecut):
        if molecule:
            Nb = 'no_Nb'
            molecule = 'molecule'
        else:
            molecule = 'suface_molecule'
        for folder in os.listdir(directory):
            path = f'{directory}/{folder}/final.e'
            e    = self.get_e(path)
            self._list.append(structure(energy      = e,
                                        name        = folder+add_name,
                                        ecut        = ecut,
                                        folder_name = path,
                                        Nb          = Nb,
                                        molecule    = molecule,
                                        functional  = func,
                                        software    = soft,
                                        pseudopot   = pp,
                                        ))
    def run(self):
        here = os.getcwd()
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/molecules',func=FUNCTIONAL.PBED3,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=True,Nb=False,ecut=400)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/molecules',func=FUNCTIONAL.PBE,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=True,Nb=False,ecut=500)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/QE/molecules/PP_samira',func=FUNCTIONAL.PBED3,soft=SOFTWARE.QE,pp=PP.SAMIRA,molecule=True,Nb=False,ecut=680)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/QE/molecules/PP_lanna',func=FUNCTIONAL.PBED3,soft=SOFTWARE.QE,pp=PP.LANNA,molecule=True,Nb=False,ecut=680)

        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/results/novdw',func=FUNCTIONAL.PBE,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=False,Nb='pv',ecut=400)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/results/vdw',func=FUNCTIONAL.PBED3,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=False,Nb='pv',ecut=400)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/results/Nb_sv_novdw',func=FUNCTIONAL.PBE,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=False,Nb='sv',ecut=400)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/results/Nb_sv_vdw',func=FUNCTIONAL.PBED3,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=False,Nb='sv',ecut=400)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/results/correct_100_vdw',func=FUNCTIONAL.PBED3,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=False,Nb='pv',add_name='_correct_100',ecut=400)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/VASP/results/500eV',func=FUNCTIONAL.PBE,soft=SOFTWARE.VASP,pp=PP.SAMIRA,molecule=False,Nb='pv',ecut=500)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/QE/PP_samira/results/novdw',func=FUNCTIONAL.PBE,soft=SOFTWARE.QE,pp=PP.SAMIRA,molecule=False,Nb='sv',ecut=680)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/QE/PP_samira/results/vdw',func=FUNCTIONAL.PBED3,soft=SOFTWARE.QE,pp=PP.SAMIRA,molecule=False,Nb='sv',ecut=680)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/QE/PP_lanna/results/novdw',func=FUNCTIONAL.PBE,soft=SOFTWARE.QE,pp=PP.LANNA,molecule=False,Nb='sv',ecut=680)
        self.iterate_folders(directory='/home/tiago.ferreiragoncal/Lanna_test/QE/PP_lanna/results/vdw',func=FUNCTIONAL.PBED3,soft=SOFTWARE.QE,pp=PP.LANNA,molecule=False,Nb='sv',ecut=680)




if __name__ == '__main__':
    properties=get_properties()
    properties.run()
    for i in properties._list:
        pass
        #print(i)


    OOH_prop   = properties.call_obj('OOH','pv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP,400)
    clean_prop = properties.call_obj('clean','pv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP,400)
    H2O_prop   = properties.call_obj('H2O','no_Nb','molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP,400)
    H2_prop    = properties.call_obj('H2','no_Nb','molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP,400)

    print(OOH_prop)
    print(clean_prop)
    print(H2O_prop)
    print(H2_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy
    H2O   = H2O_prop.energy
    H2    = H2_prop.energy

    vpvpd = (OOH-(clean+2*H2O-(3/2)*H2))
    

    OOH_prop   = properties.call_obj('OOH','pv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,400)
    clean_prop = properties.call_obj('clean','pv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,400)

    print(OOH_prop)
    print(clean_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy

    vpvp = (OOH-(clean+2*H2O-(3/2)*H2))

    OOH_prop   = properties.call_obj('OOH','sv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP,400)
    clean_prop = properties.call_obj('clean','sv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP,400)

    print(OOH_prop)
    print(clean_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy

    vsvpd = (OOH-(clean+2*H2O-(3/2)*H2))

    OOH_prop   = properties.call_obj('OOH','sv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,400)
    clean_prop = properties.call_obj('clean','sv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,400)

    print(OOH_prop)
    print(clean_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy

    vsvp = (OOH-(clean+2*H2O-(3/2)*H2))

    #OOH_prop   = properties.call_obj('OOH_correct_100','pv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP)
    #clean_prop = properties.call_obj('clean_correct_100','pv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.VASP)
    #
    #print(OOH_prop)
    #print(clean_prop)
    #
    #OOH   = OOH_prop.energy
    #clean = clean_prop.energy
    #
    #vc100pvpd = (OOH-(clean+2*H2O-(3/2)*H2))

    OOH_prop   = properties.call_obj('OOH','sv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.QE,680)
    clean_prop = properties.call_obj('clean','sv','surface_molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.QE,680)
    H2O_prop   = properties.call_obj('H2O','no_Nb','molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.QE,680)
    H2_prop    = properties.call_obj('H2','no_Nb','molecule',FUNCTIONAL.PBED3,PP.SAMIRA,SOFTWARE.QE,680)

    print(OOH_prop)
    print(clean_prop)
    print(H2O_prop)
    print(H2_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy
    H2O   = H2O_prop.energy
    H2    = H2_prop.energy

    qsvpd_s = (OOH-(clean+2*H2O-(3/2)*H2))    
   
    OOH_prop   = properties.call_obj('OOH','sv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.QE,680)
    clean_prop = properties.call_obj('clean','sv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.QE,680)

    print(OOH_prop)
    print(clean_prop)
    print('####################')
   
    OOH = OOH_prop.energy
    clean = clean_prop.energy
    
    qsvp_s = (OOH-(clean+2*H2O-(3/2)*H2))

    OOH_prop   = properties.call_obj('OOH','sv','surface_molecule',FUNCTIONAL.PBED3,PP.LANNA,SOFTWARE.QE,680)
    clean_prop = properties.call_obj('clean','sv','surface_molecule',FUNCTIONAL.PBED3,PP.LANNA,SOFTWARE.QE,680)
    H2O_prop   = properties.call_obj('H2O','no_Nb','molecule',FUNCTIONAL.PBED3,PP.LANNA,SOFTWARE.QE,680)
    H2_prop    = properties.call_obj('H2','no_Nb','molecule',FUNCTIONAL.PBED3,PP.LANNA,SOFTWARE.QE,680)

    print(OOH_prop)
    print(clean_prop)
    print(H2O_prop)
    print(H2_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy
    H2O   = H2O_prop.energy
    H2    = H2_prop.energy

    qsvpd_l = (OOH-(clean+2*H2O-(3/2)*H2))

    OOH_prop   = properties.call_obj('OOH','sv','surface_molecule',FUNCTIONAL.PBE,PP.LANNA,SOFTWARE.QE,680)
    clean_prop = properties.call_obj('clean','sv','surface_molecule',FUNCTIONAL.PBE,PP.LANNA,SOFTWARE.QE,680)

    print(OOH_prop)
    print(clean_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy

    qsvp_l = (OOH-(clean+2*H2O-(3/2)*H2))

   
    OOH_prop   = properties.call_obj('OOH','pv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,500)
    clean_prop = properties.call_obj('clean','pv','surface_molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,500)
    H2O_prop   = properties.call_obj('H2O','no_Nb','molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,500)
    H2_prop    = properties.call_obj('H2','no_Nb','molecule',FUNCTIONAL.PBE,PP.SAMIRA,SOFTWARE.VASP,500)

    print(OOH_prop)
    print(clean_prop)
    print(H2O_prop)
    print(H2_prop)
    print('####################')

    OOH   = OOH_prop.energy
    clean = clean_prop.energy
    H2O   = H2O_prop.energy
    H2    = H2_prop.energy

    vpvp_500 = (OOH-(clean+2*H2O-(3/2)*H2))

    print('\n\nReaction energies:\n')
    print(f'VASP Nb_pv PBE-D3 OOH reaction {vpvpd:.3f} eV     diff from reference:{vpvpd-vpvpd:.3f}  eV')
    print(f'VASP Nb_pv PBE OOH reaction {vpvp:.3f} eV         diff from reference:{vpvp-vpvpd:.3f}  eV')
    print(f'VASP Nb_sv PBE-D3 OOH reaction {vsvpd:.3f} eV     diff from reference:{vsvpd-vpvpd:.3f}  eV')
    print(f'VASP Nb_sv PBE OOH reaction {vsvp:.3f} eV         diff from pbe pv:{vsvp-vpvp:.3f}  eV')
    print(f'VASP Nb_sv PBE OOH reaction (500eV) {vpvp_500:.3f} eV         diff from pbe pv:{vpvp_500-vpvp:.3f}  eV')
    #print(f'VASP (correct 100) Nb_pv PBE-D3 OOH reaction {vc100pvpd:.3f} eV      diff from reference:{c100pvpd-vpvpd:.3f}  eV')  ####RUNNING####
    print('')
    print(f'QE SAMIRA PBE-D3 OOH reaction {qsvpd_s:.3f} eV     diff from reference:{qsvpd_s-qsvpd_s:.3f}  eV')
    print(f'QE SAMIRA PBE OOH reaction {qsvp_s:.3f} eV     diff from reference:{qsvp_s-qsvpd_s:.3f}  eV')
    print(f'QE LANNA PBE-D3 OOH reaction {qsvpd_l:.3f} eV     diff from reference:{qsvpd_l-qsvpd_s:.3f}  eV')
    print(f'QE LANNA PBE OOH reaction {qsvp_l:.3f} eV     diff from SAMIRA PBE:{qsvp_l-qsvp_s:.3f}  eV')
    print('')
    print(f'software difference:')
    print(f'vasp_sv QE(SAMIRA): {qsvpd_s-vsvpd:.3f} eV')
    print(f'vasp_sv QE(LANNA):  {qsvpd_l-vsvpd:.3f} eV')




