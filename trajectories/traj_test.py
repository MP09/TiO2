import ase
from ase.io.trajectory import Trajectory

traj = Trajectory('slab.traj')

A = traj[0]


for jj in A:
    print(jj)