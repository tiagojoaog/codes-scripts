#! /usr/bin/env python

import ase
from ase.io import *
import sys,os
import numpy as np
from ase.slab_tiago import slab


atoms=read("in.traj")

# bond_threshold=0.3 is the default value. If you want more structures to be generated, increase it.
slab=slab(atoms,miller=(1,2,1),repeat=(1,1,1),vacuum=18,layers=4,frozen_layers=2)

