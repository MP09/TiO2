from ase import Atoms
from ase.neighborlist import neighbor_list
import numpy as np

cell = np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]])
pos = np.array([[0, 0, 0]])

A = Atoms('C', scaled_positions=pos, cell=cell, pbc=[True, True, True])
#B = A.repeat(3)
#B.edit()

I, J, D = neighbor_list('ijd', A, 10)

print()