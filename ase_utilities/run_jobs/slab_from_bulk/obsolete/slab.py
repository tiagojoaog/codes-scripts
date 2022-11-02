#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
from ase.lattice.spacegroup import crystal
from ase.io import *
import sys
import os
import numpy as np
from ase.build import surface
from ase.build import bulk
from ase.visualize import view



bulk_traj=read("SnO2bulk.traj")
slab_traj=ase.build.surface(bulk_traj, (1,1,0),4,vacuum=8.5,periodic=True)
write('in.traj',slab_traj)

