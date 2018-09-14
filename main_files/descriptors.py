import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from itertools import permutations
from helper_functions import veclen
from scipy.special import sph_harm

def BehParCutOff(r, r_c):
    return (r<=r_c)*0.5*(1+np.cos(np.pi*r/r_c))

def BehParRadial(atoms, eta, r_c):
    num_atoms = A.get_number_of_atoms()
    num_params = len(eta)
    F = np.zeros((num_atoms, num_params))
    
    r_ij = atoms.get_all_distances()
    for ii in range(num_atoms):
        R = r_ij[ii, r_ij[ii, :] > 0].reshape(num_atoms-1, 1)
        F[ii, :] = np.sum(np.exp(-eta*R**2/r_c**2)*BehParCutOff(R, r_c), axis=0)
    return F

def BehParAngular(atoms, eta, xi, r_c, lambd=[1, -1]):
    num_atoms = A.get_number_of_atoms()
    num_params = 2*(len(eta)*len(xi))
    F = np.zeros((num_atoms, num_params))
    
    fc = BehParCutOff

    r_ij = atoms.get_all_distances()
    for ii in range(num_atoms):
        for jj in range(num_atoms):
            for kk in range(num_atoms):
                if kk != ii and jj != kk and ii != jj:
                    theta = atoms.get_angle(jj, ii, kk)*np.pi/180
                    c = 0
                    fac = fc(r_ij[ii, jj], r_c)*fc(r_ij[ii, kk], r_c)*fc(r_ij[jj, kk], r_c)
                    for e in eta:
                        for x in xi:
                            for lamb in lambd:
                                F[ii, c] += 2**(1-x)*(1-lamb*np.cos(theta))**x*np.exp(-e*(r_ij[ii, jj]**2+r_ij[ii, kk]**2+r_ij[jj, kk]**2)/r_c**2)*fac
                                c += 1
    return F                    

def BehPar(atoms, eta, xi, rc):
#    print(eta)
#    print(xi)
#    print(rc)
    I, J, dists, D = neighbor_list('ijdD', atoms, rc)
    dists = dists[:, np.newaxis]

    num_radial = len(eta)
    num_angular = len(xi)*2
    num_atoms = len(atoms)
    F = np.zeros((num_atoms, num_radial+num_angular))
    
    # Adjust angular parameters:
    lamb = np.zeros((num_angular))
    Xi = np.zeros((num_angular))
    c = 0
    for x in xi:
        for ii in [-1, 1]:
            Xi[c] = x
            lamb[c] = ii
            c += 1

    # Radial functions given by:
    for i, j, d in zip(I, J, dists):
        F[i, 0:num_radial] += np.exp(-eta*d**2/rc**2)*BehParCutOff(d, rc)
    
    # Angular functions:
    eta_ang = 0.005
    for i in range(num_atoms):
        neigh_mask = I==i
        for Rij, rij, j in zip(D[neigh_mask], dists[neigh_mask].ravel(), J[neigh_mask]):
            for Rik, rik, k in zip(D[neigh_mask], dists[neigh_mask].ravel(), J[neigh_mask]):
                if j < k:
                    Rjk = Rik-Rij
                    rjk = np.sqrt(Rjk@Rjk)
                    theta = (Rij@Rik)/(np.sqrt(Rij@Rij)*np.sqrt(Rik@Rik))
                    F[i, num_radial::] += (1+lamb*theta)**Xi*np.exp(-eta_ang*(rij**2+rik**2+rjk**2)/rc**2)*BehParCutOff(rij, rc)*BehParCutOff(rik, rc)*BehParCutOff(rjk, rc)
    F[:, num_radial::] *= 2**(1-Xi)
    
    # Atomic number as last entry:
    Fz = atoms.get_atomic_numbers().reshape(num_atoms, 1)
    
    
    #F[: -1] = [atom.get_atomic_number() for atom in atoms]

    F = np.concatenate((Fz, F), axis=1)

    return F

def SivaDescriptor(atoms, eta, xi, rc, L=2):
    """
    Calculates the descriptors of an atoms object according to 
    'S. Jindal, S. Chiriki, Spherical Harmonics based descriptor..'
    -- atoms: Atoms objects
    -- eta: Parameters for radial functions
    -- xi: Parameters for angular functions
    -- rc: Cutoff
    -- L: Maximum l used in expansion of spherical harmonics.

    Returns local feature vectors of each atom in the structure:
    1st entry: Chemical species.
    Middle: Radial functions
    Last: Angular functions
    """


    # Obtain data from atoms object
    I, J, dists, Dvecs = neighbor_list('ijdD', atoms, rc)
    num_atoms = atoms.get_number_of_atoms()
    num_eta = len(eta)
    num_xi = len(xi)
    num_l = L+1; Ls = [l for l in range(num_l)]
    num_m = 2*L+1

    # Setup storage arrays:
    Fr = np.zeros((num_atoms, num_eta)) # Radial
    Fa = np.zeros((atoms.get_number_of_atoms(), num_xi, num_l)) # Angular

    # Calculate radial functions:
    for i, d in zip(I, dists):
        Fr[i, :] += np.exp(-eta*d**2)
    
    # Calculate angular functions:
    c = np.zeros((num_atoms, num_xi, num_l, num_m), dtype=np.complex)
    for i, d, r in zip(I, dists, Dvecs):

        pol = np.arccos(r[2]/d)    # Polar angle
        azi = np.arctan(r[1]/r[0]) # Azimuthal angle
        fc = BehParCutOff(d, rc)

        for xidx, x in enumerate(xi):
            ans = np.exp(-x*d**2)*fc
            for lidx, l in enumerate(range(L+1)):
                for midx, m in enumerate(range(-l, l+1)):        
                    c[i, xidx, lidx, midx] = ans*sph_harm(m, l, azi, pol)

    for i in range(num_atoms):
        for lidx in range(num_l):
            for xidx in range(num_xi):
                Fa[i, xidx, lidx] += 4*np.pi/(2*Ls[lidx]+1)*(np.vdot(c[i, xidx, lidx, :],c[i, xidx, lidx, :])).real
    
    Fa = Fa.reshape(num_atoms, num_xi*num_l)
    Z = atoms.get_atomic_numbers()[:, np.newaxis]
    return np.append(Z, np.append(Fr, Fa, axis=1), axis=1)







