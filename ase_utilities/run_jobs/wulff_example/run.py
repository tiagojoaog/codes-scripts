#! /usr/bin/env python3

from wulffpack import (SingleCrystal,
                       Decahedron,
                       Icosahedron)
from ase.build import bulk
from ase.io import read,write


surface_energies = {(1, 0, 0): 1.0,
                    (1, 1, 1): 0.9,
                    (1, 1, 0): 1.1}

prim = read('LaAlO3_mp-5304_computed.cif')
particle = SingleCrystal(surface_energies,
                      primitive_structure=prim,
                      natoms=600)
#particle.view()
write('nanoparticle.xyz', particle.atoms)
