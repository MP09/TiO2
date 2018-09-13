import numpy as np
from timeit import default_timer as dt

def fireworks(struct, idx, R):
    """
    Moves the positions of the atoms given by idx according to the fireworks pertubation
    -- struct: Structure object
    -- idx: Indicies of atoms to move.
    -- R: Radius of sphere to move from center of mass.
    """
    atoms = struct.atoms
    com = atoms.get_center_of_mass()

    t0 = dt()
    for ii in idx:
        not_found = True
        while not_found:
            xyz = UnitSphereDist(R=R)
            xyz += com
            if (xyz > struct.min_dims).all():
                if (xyz < struct.max_dims).all():
                    if ProximityCheck(atoms, xyz, ii):
                        atoms.positions[ii, :] = xyz
                        not_found = False

def random_center(struct, idx, R):
    """
    Moves the atoms given by idx to a random position above another atom.
    -- struct: Structure object
    -- idx: Indicies of atoms to move.
    -- R: Radius of sphere to move from center of mass.
    """
    atoms = struct.atoms
    num_atoms = atoms.get_number_of_atoms()
    allowed_idx = [ii for ii in range(num_atoms) if ii not in idx]

    zrange = R #(struct.max_dims[2]-struct.min_dims[2])
    for ii in idx:
        not_found = True
        while not_found:
            jj = np.random.choice(allowed_idx)
            xyz = atoms[jj].position.copy()
            xyz[2] += np.random.rand()*zrange
            if (xyz > struct.min_dims).all():
                if (xyz < struct.max_dims).all():
                    if ProximityCheck(atoms, xyz, ii):
                        atoms.positions[ii, :] = xyz
                        not_found = False

def random_position(struct, idx, *args):
    """
    Moves the atoms given by idx to random positions within the constraints set by the min/max dimensions of the structure.
    """
    atoms = struct.atoms

    bias = np.array([0, 0, 5])
    for ii in idx:
        not_found = True
        while not_found:
            xyz = struct.min_dims + np.random.rand(3)*(struct.max_dims-struct.min_dims-bias)
            if ProximityCheck(atoms, xyz, ii):
                atoms.positions[ii, :] = xyz
                not_found = False

# Helper functions
def UnitSphereDist(R=1):
    xyz = np.random.uniform(low=-1, high=1, size=3); xyz *= 1/np.linalg.norm(xyz)*R*np.random.rand(1)
    return xyz

def ProximityCheck(atoms, xyz, idx, min_dist=0.35):
    dists = {'O':0.66, 'Ti': 1.60}
    d1 = dists[atoms[idx].symbol]
    for jj in range(atoms.get_number_of_atoms()):
        if jj != idx:
            d2 = dists[atoms[jj].symbol]
            if np.linalg.norm(atoms[jj].position-xyz) < 0.7*(d1+d2):
                #print(jj, 'kage')
                return False
    return True

def UnitCellCheck(atoms):
    """
    Check that all atoms are within the unit cell, assumes a cubic cell.
    """
    uc = atoms.get_cell()
    max_dims = np.diag(uc)
    for atom in atoms:
        if (atom.position > max_dims).any():
            return False
    return True

pertubation_dict = {0:['fireworks', 'FW', fireworks], 1:['random center', 'RC', random_center], 2:['random position', 'RP', random_position]}
















