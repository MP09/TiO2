from main import TiO2, autobag
import numpy as np
from descriptors import SivaDescriptor

a = TiO2(local=True)

a.add_atoms(num_units=1)


eta = np.array([0, 1])
xi = np.array([])
rc = 10

params = [eta, xi, rc, 1]

A = autobag(params)

A.add_structure(a.atoms, 1)
A.add_structure(a.atoms, 1)
A.init_kmeans()
A.calc_global_features()
A.calc_local_energy()
print(A.pick_atom(2))
