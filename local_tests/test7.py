import numpy as np
from ase.io import Trajectory
from basin_hoping import basin_hopping, TiO2
from pertubations import fireworks, random_center, random_position

A = TiO2(local=True)
A.add_atoms(num_units=5)
atoms = A.atoms

num_atoms = len(atoms)

traj = Trajectory('test7', mode='w')
for jj in range(100):

    idx = [np.random.randint(low=0, high=num_atoms)]
    random_position(A, idx)

    traj.write(A.atoms)




