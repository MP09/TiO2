import numpy as np

from basin_hoping import basin_hopping, TiO2
from pertubations import fireworks, UnitCellCheck

A = TiO2(local=True)
A.add_atoms(num_units=1)

atoms = A.atoms
#atoms.edit()
#atoms.edit()

UnitCellCheck(atoms)




#Ba = basin_hopping()

#Ba.start()