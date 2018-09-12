import numpy as np

from basin_hoping import basin_hopping, TiO2
from pertubations import fireworks

A = TiO2(local=True)
A.add_atoms(num_units=1)

atoms = A.atoms
atoms.edit()

fireworks(A, [12, 13, 14], 10)
atoms.edit()




#Ba = basin_hopping()

#Ba.start()