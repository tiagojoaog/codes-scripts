#!/usr/bin/env python3

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List
from matplotlib import pyplot as plt
from ase.io import read

import sys, os, warnings
import numpy as np

'''

Example script to compute surface energies
for a surface that has the same facet but 
different terminations.

'''


class Functional(Enum):
    """Functional chosen."""

    PBE         = auto()
    PBE_SPINPOL = auto()
    RPBE_SP     = auto()
    RPBE_OPT    = auto()


class NoEnergyOrVibsError(Exception):
    """ Custom error that is raised if no energy or vibrations were extracted
        successfully."""

    def __init__(self, energy: float, vibs: float, structure: str,
            functional: Functional, location: str, warn: bool = False) -> None:

        self.energy     = energy
        self.vibs       = vibs
        self.structure  = structure
        self.functional = functional
        self.location   = location
        self.warn       = warn
    
        message         = f"""Could not extract energies and/or vibs. 

                              Structure : {self.structure}
                              Folder    : {self.location}
                              Functional: {self.functional}
                              Energy    : {self.energy}
                              Vibrations: {self.vibs}"""
        if self.warn:     
            message        += """
                              Note      : Total energy and/or vibration correction will be considered 0.
                                          A placeholder can be defined later.
                          """
        self.message    = message
        super().__init__(message)

class NoEnergyOrVibsWarning(UserWarning):
        """ Custom Warning that is raised if no energy or vibrations were extracted
            successfully."""

        pass


@dataclass(order=True) #frozen=False by default, if database is read only (not allowed to change values assigned), use frozen=True
class Structure:
    """ Model of the catalyst (surface, zeolite, molecule) with an adsorbate
        along with the functional chosen for total energy calculation and 
        vibrational analysis using Rigid rotator/translator and harmonic
        oscillator approximation. """

    sort_index  : int   = field(init=False, repr=False)
    name        : str
    index       : str
    functional  : Functional
    location    : str
    energy      : float = None  # in eV
    vibs        : float = None  # in eV

    def __post_init__(self) -> None:
        #self.sort_index = self.total_energy   # does not work if dataclass is frozen
        object.__setattr__(self, 'sort_index', self.total_energy)  # also works if dataclass is frozen

    def __str__(self) -> str:
        return f""" 
                    {self.name}:
                      * Folder                      --> {self.location}
                      * Functional                  --> {self.functional}
                      * Electronic Energy (eV)      --> {self.energy}
                      * Vibrational Correction (eV) --> {self.vibs}
                """

    @property
    def total_energy(self) -> float:
        return self._total_energy

    @total_energy.getter
    def total_energy(self) -> float:
        """ Returns energy + vibs for a given functional.
            If energy, or vibs is None, that value will
            be considered zero. """

        _e = self.energy
        _v = self.vibs

        if _e == None:
            _e = 0.0
        if _v == None:
            _v = 0.0

        return _e + _v

@dataclass
class Surface:
    """ Surface type definition. """

    index           : str
    functional      : Functional
    surface_area    : float      # in Angstroem
    number_of_atoms : int

    def __str__(self) -> str:
        return f" Surface {self.index} with surface area of {self.surface_area} Angstroem using {self.functional} with {self.number_of_atoms} atoms."

@dataclass
class Bulk:
    """ Bulk type definition. """

    functional      : Functional
    number_of_atoms : int
    energy          : float      # in eV

    def __str__(self) -> str:
        return f" Bulk with {self.number_of_atoms} atoms, using {self.functional}. Energy is {self.energy} eV."

def placeholder(structure: object, value : float, which : str,verbose: int) -> None:
    """ Add a placeholder for unfinished calculations. """

    if which == 'energy':
        structure.energy = value
    elif which == 'vibs':
        structure.vibs = value
    
    if verbose > 1:
        print(f""" 
                        A placeholder was added to {which} values.
                          * Structure  --> {structure.name} 
                          * Functional --> {structure.functional}
                          * Folder     --> {structure.location}
              """)

def extract_energy(name : str) -> float:
    """ Extracts electronic energies from files."""

    if os.path.exists(name):
        with open(name) as val:
            num = [i for i in val]
            num = float(num[0])
    else:
        num = None
    return num
            
def get_surfaces(L: List, location: str) -> None:
    """ Obtains surface information. """

    folders  = os.listdir(location)
  
    to_remove= []
    for i in folders: # Remove folders that are not folders from list and rpbe_opt
        if len(i.split("_")) == 1 or i == 'rpbe_opt':
            to_remove.append(i)
    for i in to_remove:
        for j in folders:
            if i == j:
                folders.remove(i)

    folders+=['rpbe_opt/'+i for i in os.listdir(location+'rpbe_opt')]
    for string in folders:   #read surface folder, extract natoms, surface area, index
        if string.split("_")[0] == 'surface':
            functional = Functional.PBE
        elif string.split("_")[0] == 'rpbe':
            functional = Functional.RPBE_OPT

        try:
            os.chdir(location+string)
            atoms        = read('slab_original/in.traj')
            natoms       = len(atoms.get_positions())
            index        = string.split("_")[-1]
            surface_area = np.linalg.det([[atoms.get_cell()[0][0],atoms.get_cell()[0][1]],[atoms.get_cell()[1][0],atoms.get_cell()[1][1]]])
        
            L.append(Surface(index=index,functional=functional,surface_area=surface_area,number_of_atoms=natoms))
        except FileNotFoundError:
            pass #for now, will remove after I have the folders

def get_structures(L: List, surfaces: List, structures: List, location: str) -> None:
    """ Obtains structure information. """

    pbe_rpbe = ['.','rpbe_opt']

    for PorR in pbe_rpbe:
        for surf in surfaces:
            subfolder = f'{PorR}/surface_{surf.index}'
            os.chdir(location+subfolder)
            calculations = os.listdir(location+subfolder)
            for calculation in calculations:
                folder = calculation.split("_")
                if len(folder) > 1:
                    if folder[0] == 'slab' and PorR == '.' and not surf.functional == Functional.RPBE_OPT:
                        final_location = subfolder+'/'+calculation
                        os.chdir(location+final_location)

                        pbe         = extract_energy("final.e")
                        pbe_spinpol = extract_energy("spin_pol/final.e")
                        #rpbe_sp     = extract_energy("rpbe/final.e")

                        structures.append(Structure(name=calculation,index=surf.index, functional=Functional.PBE, location=location+final_location, energy=pbe, vibs=0))
                        structures.append(Structure(name=calculation,index=surf.index, functional=Functional.PBE_SPINPOL, location=location+final_location, energy=pbe_spinpol, vibs=0))
                        #structures.append(Structure(name=calculation,index=surf.index, functional=Functional.RPBE_SP, location=location+final_location, energy=pbe_spinpol, vibs=0))

                    elif folder[1] == 'adsorbed' and PorR == '.' and not surf.functional == Functional.RPBE_OPT:
                        pass #to change at a later stage
                    
                    elif PorR == 'rpbe_opt' and folder[0] == 'slab' and surf.functional == Functional.RPBE_OPT:
                        #rpbe_opt   = extract_energy('final.e')
                        #structures.append(Structure(name=calculation,index=surf.index, functional=Functional.RPBE_OPT, location=location+final_location, energy=rpbe_opt, vibs=0))
                        pass

                    elif folder[1] == 'adsorbed' and PorR == 'rpbe_opt' and surf.functional == Functional.RPBE_OPT:
                        pass #change this as well


def RaiseErrorOrWarning(no_placeholders: bool, verbose: str, Structures: object) -> None:
    """ Plots surface energy in meV/Angstroem. """

    for i in Structures:
        if (i.energy == None or i.vibs == None) and no_placeholders:
            raise NoEnergyOrVibsError(energy   = i.energy,
                                    vibs       = i.vibs,
                                    structure  = i.name,
                                    functional = i.functional,
                                    location   = i.location)
        elif (i.energy == None or i.vibs == None) and verbose>2:
            warnings.warn(f"""Could not extract energies and/or vibs.
                            Structure : {i.name}
                            Folder    : {i.location}
                            Functional: {i.functional}
                            Energy    : {i.energy}
                            Vibrations: {i.vibs}""",NoEnergyOrVibsWarning)

    

def plot_surface_energy(index: str, Structures: object, Surfaces: object, Bulk: List, functional: object) -> List:
    """ Plots surface energy in meV/Angstroem. """

    plot_y    = []
    labels    = []
    bulk_PBE  = Bulk[0]
    bulk_RPBE = Bulk[1]

    for i in Surfaces:
        if functional == Functional.PBE_SPINPOL or functional == Functional.RPBE_SP:
            surf_functional=Functional.PBE
        else:
            surf_functional=functional
        for j in Structures:
            if (j.functional == functional and i.functional == surf_functional) and j.index == i.index and i.index == index:
                plot_y.append(((j.energy - bulk_PBE.energy*(i.number_of_atoms/bulk_PBE.number_of_atoms))*1000/(2*i.surface_area))) #displays in meV/Angstrom
                labels.append(j.name)

    plot_x = [i for i in range(len(plot_y))]

    if functional.name == 'PBE':
        func = 'PBE-D3'
    elif functional.name == 'PBE_SPINPOL':
        func= 'PBE-D3_SPINPOL'
    else:
        func = functional.name

    plt.plot(plot_x,plot_y,label=func, marker='.')
    
    return labels, plot_x

def main() -> None:
    """ Main function. """

    no_placeholders = False
    verbose     = 3
    bulk_PBE    = Bulk(functional=Functional.PBE, number_of_atoms=24, energy=-181.05002053)
    bulk_RPBE   = Bulk(functional=Functional.RPBE_OPT, number_of_atoms=24, energy=-178.43966238)
    bulk        = [bulk_PBE, bulk_RPBE]
    location    = "/home/tiago.ferreiragoncal/H2O2_project/BiVO4/surfaces/conventional/"
    here        = os.getcwd()
    Structures  = []
    Surfaces    = []

    get_surfaces(L=Surfaces, location=location)

    for index, value in enumerate(sys.argv):

        if value == "-strict":
            no_placeholders = True  #Used to check incomplete calculations

        elif value == "-v":
            verbose = int(sys.argv[index+1]) # v=1 simple print v=2 intermediate print (shows where placeholders have been used) v>=3 full print (also shows warnings)

    get_structures(L=Structures, surfaces=Surfaces, structures=Structures, location=location)
    RaiseErrorOrWarning(no_placeholders=no_placeholders, verbose=verbose, Structures=Structures) 
    os.chdir(here)

    labels,x=plot_surface_energy(index='020', Structures=Structures, Surfaces=Surfaces, Bulk=bulk, functional=Functional.PBE)
    labels,x=plot_surface_energy(index='020', Structures=Structures, Surfaces=Surfaces, Bulk=bulk, functional=Functional.PBE_SPINPOL)
    plt.title('BiVO4 020 surface')
    #plt.xlabel('Structure')
    plt.ylabel(r'Surface energy (meV/$\AA$)')
    plt.legend(loc='best')
    plt.xticks(x, labels, rotation='vertical')
    plt.tight_layout()
    plt.savefig('surf_020.pdf')
    plt.close()

    fig, ax1 = plt.subplots()
    fig.set_size_inches(10.0, 5)
    labels,x=plot_surface_energy(index='121', Structures=Structures, Surfaces=Surfaces, Bulk=bulk, functional=Functional.PBE)
    labels,x=plot_surface_energy(index='121', Structures=Structures, Surfaces=Surfaces, Bulk=bulk, functional=Functional.PBE_SPINPOL)
    plt.title('BiVO4 121 surface')
    #plt.xlabel('Structure')
    plt.ylabel(r'Surface energy (meV/$\AA$)')
    plt.legend(loc='best')
    plt.xticks(x, labels, rotation=90)
    plt.tight_layout()
    plt.savefig('surf_121.pdf')
    plt.close()

    labels,x=plot_surface_energy(index='211', Structures=Structures, Surfaces=Surfaces, Bulk=bulk, functional=Functional.PBE)
    labels,x=plot_surface_energy(index='211', Structures=Structures, Surfaces=Surfaces, Bulk=bulk, functional=Functional.PBE_SPINPOL)
    plt.title('BiVO4 211 surface')
    #plt.xlabel('Structure')
    plt.ylabel(r'Surface energy (meV/$\AA$)')
    plt.legend(loc='best')
    plt.xticks(x, labels, rotation=90)
    plt.tight_layout()
    plt.savefig('surf_211.pdf')
    plt.close()


if __name__ == "__main__":
    main()
    sys.exit(0)


