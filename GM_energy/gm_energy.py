from basin_hoping import TiO2
from ase.io import Trajectory, write, read

A = TiO2(local=False)

A.forcemax = 0.005

#A.add_atoms()

gm_atoms = Trajectory(A.prefix+'/trajectories/gm.traj')[0]

A.set_atoms(gm_atoms)

uc = A.atoms.get_cell()

#A.atoms[12].position = A.atoms[12].position + 10*uc[0, :]

#A.atoms.edit()

E = A.relax_structure_V2()


#E = A.get_energy()
print(E)
A.atoms.edit()


