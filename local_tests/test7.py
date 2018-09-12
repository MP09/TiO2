import numpy as np

from basin_hoping import basin_hopping, TiO2
from pertubations import fireworks, random_center

A = TiO2(local=True)
A.add_atoms(num_units=5)
atoms = A.atoms
#atoms.edit()

random_center(A, [12], 15)
atoms.edit()

