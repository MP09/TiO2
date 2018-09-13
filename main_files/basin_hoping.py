#from ase import MPI
from ase.calculators.dftb import Dftb
from ase.optimize import BFGS, QuasiNewton
from ase.optimize.sciopt import SciPyFminCG
from ase import Atom
from ase.io import Trajectory, write, read

from sklearn.cluster import KMeans

import numpy as np

from descriptors import SivaDescriptor, BehPar
from helper_functions import NormGeoSeries
from pertubations import fireworks, UnitCellCheck, random_center

from os.path import expanduser

from timeit import default_timer as dt
from ase.ga.relax_attaches import VariansBreak

class basin_hopping:

    def __init__(self, params_structure=[], use_ML=False, traj_name='BH.traj'):
        """
        Takes a structure object and performs a Basin-hopping search for the global minimum.
        """

        if expanduser('~') == '/home/mp':
            self.local = True
        else:
            self.local = False

        # Set the structure class:
        self.structure_class = TiO2
        self.init_structure = self.structure_class(*params_structure, local=self.local)
        self.init_structure.add_atoms()

        # Information about the structure:
        self.num_atoms = self.init_structure.num_atoms
        self.num_atoms_perturb = self.num_atoms-len(self.init_structure.atoms.constraints[0].index)

        # Settings:
        self.max_iterations = 10
        self.GMenergy = -1933.30306 # Checked 20/8
        self.perturb_params = 10
        self.file_prefix = 'TiO2_{}'

        # Other stuff.
        self.W = np.cumsum(NormGeoSeries(np.arange(self.num_atoms_perturb)+1, self.num_atoms_perturb))
        
        # Autobag stuff:
        self.use_ML = use_ML
        if self.use_ML:
            eta = np.array([0.005, 0.015, 0.0230, 0.038])
            xi = np.array([0.0028, 0.0048])
            rc = 10
            L = 3
            params = [eta, xi, rc, L]
            self.autobag = autobag(params)

        # File names
        if self.use_ML:
            self.trajectory_file = 'ML_'+traj_name
        else:
            self.trajectory_file = traj_name
        self.trajectory_obj = Trajectory(self.trajectory_file, mode='w')
        self.trajectory_obj_perturb = Trajectory('P_'+self.trajectory_file, mode='w')

        self.perturb = fireworks


        # Intial printing:
        print('Settings:')
        print('Using machine learning: ' + str(self.use_ML))

    def set_perturbation(self, label):
        if label == 'fireworks':
            print('Using the ' + label + 'pertubation')
            self.perturb = fireworks
        elif label == 'random_center':
            print('Using the ' + label + 'pertubation')
            self.perturb = random_center
        else:
            print('This is bad. Check pertubation label')

    def start(self):
        """
        Starts basin hoping algorithm:
        """
     
        print('Basinhoping started', flush=True)

        # Initialize the first structure:
        cell_check = UnitCellCheck(self.init_structure.atoms)
        

        best_E = self.init_structure.relax_structure_V2()
        best = self.init_structure 

        cell_check = UnitCellCheck(best.atoms)


        # Add to autobag:
        if self.use_ML:
            self.autobag.add_structure(self.init_structure.atoms, best_E)

        # Arrays/Settings:
        niter = 0
        space = ' '*8

        while niter < self.max_iterations:
            # Make new structure:
            cand = self.structure_class(local=self.local, post_fix='{}'.format(niter))
            cand.set_atoms(best.atoms.copy())
            # Perturb new structure: 
            N = self.calc_number_of_changes()
            idx = self.pick_func(N)
            self.perturb(cand, idx, self.perturb_params)
            self.trajectory_obj_perturb.write(cand.atoms)

            # Relax the now perturbed structure:
            t0 = dt()
            cand_E = cand.relax_structure_V2(label=self.file_prefix.format(niter))
            relax_time = dt()-t0
            
            
            # Check whether a better structure was found:
            if cand_E < best_E:
                best = cand
                best_E = cand_E
                self.trajectory_obj.write(best.atoms)

            print('#'+20*'='+'Iteration {}'.format(niter)+20*'='+'#', flush=True)
            print(space+'Number of pertubed atoms: {}'.format(N))
            print(space+'Relaxation took {:6.2f} s'.format(relax_time))
            print(space+'Best energy: {:6.2f} eV'.format(best_E), flush=True)
            print(space+'Cand. energy: {:6.2f} eV'.format(cand_E), flush=True)

            # Add to autobag if using ML:
            if self.use_ML:
                self.autobag.add_structure(cand.atoms, cand_E)

            if best_E <= self.GMenergy:
                return True, niter, best_E
                print('Calculation converged!')
            niter += 1


        return False, niter, best_E
        print('Warning: Not converged')


    #def perturb(self, cand, idx):
    #    random_center(cand, idx, self.FireworkRadius)
    #    #fireworks(cand, idx, self.FireworkRadius)

    def pick_func(self, N):
        if self.use_ML:
            return self.autobag.pick_atom(N)   
        else:
            cons = self.init_structure.atoms.constraints[0].index
            pos_idx = np.array([idx for idx in range(self.num_atoms) if idx not in cons])
            return np.random.choice(pos_idx, size=N)

    def calc_number_of_changes(self):
        """
        Calculate the number of atoms to perturb according to geometric series
        """
        return np.sum(self.W<np.random.rand())+1
        
class autobag:
    def __init__(self, params):
        """
        Only works if the descriptor function takes arguments in the order: descriptor(atoms, params[0], params[1])
        -- params: List of parameters for descriptor calculation, first entry --> first argument.
        
        Currently all added structures are considered training structures.
        """
        self.params = params

        # Atom objects:
        self.structures = [] 
        self.num_structures = 0

        # Settings:
        self.descriptor_func = SivaDescriptor # Set the descriptor function:
        self.num_clusters = 10 # Number of clusters.
        self.lamb = 0.01 # Lambda factor for regression.
        
    def add_structure(self, atoms, E):
        """
        Add a structure to the autobag:
        -- atoms: ASE Atoms object.
        -- E: Energy of structure.
        """

        self.structures.append(atoms)
        self.num_structures += 1
        try:
            self.LFVs = np.append(self.LFVs, self.descriptor_func(atoms, *self.params), axis=0)
            self.E = np.append(self.E, E)
        except AttributeError:
            self.LFVs = self.descriptor_func(atoms, *self.params)
            self.num_atoms = atoms.get_number_of_atoms()
            self.E = np.array([E])

    def init_kmeans(self):
        self.KMeans = KMeans(n_clusters=self.num_clusters)

    def calc_global_features(self):
        """
        Clusters and calculates global feature vectors of all added structures.
        """
        self.S = self.KMeans.fit_predict(self.LFVs).reshape(self.num_structures, self.num_atoms)
        self.GFVs = np.zeros((self.num_structures, self.num_clusters))
        for idx, s in enumerate(self.S):
            idxs, counts = np.unique(s, return_counts=True)
            self.GFVs[idx, idxs] = counts
        
    def calc_local_energy(self):
        """
        Calculate the local energies corresponding to each cluster.
        """
        X = self.GFVs; XT = np.transpose(X)
        self.eps = np.dot(np.linalg.inv(np.dot(XT, X)+self.lamb*np.eye(self.num_clusters)),np.dot(XT, self.E))

    def pick_atom(self, N):
        """
        Pick N atoms taking into account local energy and constraints of the most recently added structure.

        Note: Only works if FixAtoms is the only constraint.
        """
        self.init_kmeans()
        self.calc_global_features()
        self.calc_local_energy()

        # Determine indicies of the atoms that can possibly be moved:
        atom_idx = []
        cons = self.structures[-1].constraints[0].index
        for n in range(N):
            pos_idx = np.array([idx for idx in range(self.num_atoms) if idx not in cons])
            num_pos = len(pos_idx) # Number of possible choices

            # Find the local energies of these atoms:
            AEs = self.eps[self.S[-1, pos_idx]]

            # Sort according to energy:
            sidx = np.argsort(AEs).astype(int); 
            AEs = AEs[sidx]; pos_idx = pos_idx[sidx]
            
            # Calculate the chance/weight of picking each:
            W = np.cumsum(NormGeoSeries(np.arange(num_pos)+1, num_pos))
            idx = pos_idx[num_pos-np.sum(W>np.random.rand())]
            cons = np.append(cons, idx); atom_idx.append(idx)
        return atom_idx

class TiO2:

    def __init__(self, calc_type='Hotbit', local=False, post_fix=''):

        # Filepaths:
        self.local = local
        if self.local == False:
            self.prefix = '/home/machri/TiO2'
        else:
            self.prefix = '/home/mp/Dropbox/Python/TiO2'

        self.postfix = post_fix
        self.log_file = None#'TiO2_'+self.postfix+'.log'
        self.trajectory_file = None #'TiO2'+self.postfix+'.traj'
    
        # Settings:
        self.z_bias = 3 # Minimum distance of extra atoms to cell boundary in z.
        self.forcemax = 0.075 # Maximim force during local relaxation

        # ASE Atoms object:
        self.atoms = Trajectory(self.prefix+'/trajectories/slab.traj')[0]
        self.num_atoms_slab = len(self.atoms)
        self.min_z = max(self.atoms.get_positions()[:, 2])
        self.max_x = self.atoms.get_cell()[0, 0]
        self.max_y = self.atoms.get_cell()[1, 1]
        self.max_z = self.atoms.get_cell()[2, 2] #-self.min_z-self.z_bias
        
        self.min_dims = np.array([0, 0, self.min_z])
        self.max_dims = np.array([self.max_x, self.max_y, self.max_z])


        # Specifies which calculator to use: Hotbit or DFTB+
        self.calc_type = calc_type

    def set_atoms(self, atoms):
        """
        Sets the self.atoms object to the given atoms object.
        """
        self.atoms = atoms

    def add_atoms(self, num_units=5):
        """
        Add extra atoms to the surface.
        """
    
        dims = np.array([self.max_x, self.max_y, self.max_z-self.min_z-self.z_bias])
        for jj in range(num_units):
            for atom_str in ['Ti', 'O', 'O']:
                self.atoms.append(Atom(atom_str))
                self.atoms[-1].position = np.random.rand(3)*dims
                self.atoms[-1].position[2] += self.min_z
        self.num_atoms = self.atoms.get_number_of_atoms()

    def set_calculator(self, label):
        self.atoms.set_calculator(Dftb(label, atoms=self.atoms, kpts=(2, 1, 1),
                    Hamiltonian_Eigensolver=r'RelativelyRobust {}',    
                    Hamiltonian_MaxAngularMomentum_O='p',
                    Hamiltonian_MaxAngularMomentum_Ti='d',
                    Hamiltonian_Periodic='Yes',
                    Hamiltonian_SCC ='No',
                    Hamiltonian_Filling ='Fermi {',
                    Hamiltonian_Filling_empty= 'Temperature [Kelvin] = 0.000000'))

        #tables = {'TiTi':'/home/machri/software/DFTB+_SK_files/matsci-0-3/Ti-Ti.skf', 
        #          'TiO':'/home/machri/software/DFTB+_SK_files/matsci-0-3/Ti-O.skf',
        #          'OO':'/home/machri/software/DFTB+_SK_files/matsci-0-3/O-O.skf'}

        #from hotbit import Hotbit
        #calc = Hotbit(tables=tables, SCC=False, width=0.001, kpts=(2, 1, 1), txt='log.txt')
        #print(dir(calc))
        #print(calc.tables)

        #from ase.calculators.emt import EMT
        #self.atoms.set_calculator(EMT())
        #print('Calculator set')
        #from gpaw import GPAW
        #calc = GPAW(mode='lcao', basis='sz', h=0.2)
        
        #self.atoms.set_calculator(calc)


    def relax_structure_DFTB(self, label='TiO2'):
        calc = Dftb(label, atoms=self.atoms,
                    run_manyDftb_steps=True,
                    Driver_='ConjugateGradient',
                    Driver_MaxForceComponent='0.05',
                    Driver_MaxSteps=1000,
                    Driver_Constraints= self.constraint_to_string(),
                    Hamiltonian_KPointsAndWeights='''SupercellFolding { 
                        2  0  0
                        0  2  0 
                        0  0  1
                        0.5  0.5  0
                        }''',
                    Hamiltonian_MaxAngularMomentum_O='p',
                    Hamiltonian_MaxAngularMomentum_Ti='d',
                    Hamiltonian_Periodic='Yes',
                    Hamiltonian_SCC ='No',
                    Hamiltonian_Filling ='Fermi {',
                    Hamiltonian_Filling_empty= 'Temperature [Kelvin] = 0.000000')

        self.atoms.set_calculator(calc)
        t0 = dt()
        calc.calculate(self.atoms)
        print('Relaxation took: {} s'.format(dt()-t0))
        self.atoms = read('geo_end.gen')
        self.atoms.edit()

    def relax_structure(self, label='TiO2'):
        """
        Relax structure using DFTB. 
        """
        self.set_calculator(label)

        #print('Relaxation starting..', flush=True)
        dyn = BFGS(self.atoms, trajectory=self.trajectory_file, logfile=self.log_file)    
        dyn.run(fmax=self.forcemax)
        self.energy = self.atoms.get_potential_energy()
        #print('Relaxation finished', flush=True)

        #self.atoms.edit()
        #final = read('geo_end.gen')
        #write(self.trajectory_file, final)

    def relax_structure_V2(self, label='TiO2'):
        self.set_calculator(label)
        niter = 0

        while (self.atoms.get_forces()**2).sum(axis = 1).max()**0.5 > self.forcemax and niter < 10:                             
            #dyn = BFGS(self.atoms, lofile=None)                                                                                                                                                    
            #vb = VariansBreak(self.atoms, dyn, min_stdev = 0.01, N = 15)                                                   
            #dyn.attach(vb)                                                                                                
            #dyn = SciPyFminCG(self.atoms)
            dyn = QuasiNewton(self.atoms, logfile=None)
            dyn.run(fmax=self.forcemax)                                                                                      
            niter += 1

        E = self.atoms.get_potential_energy()
        self.atoms.wrap()
        return E

    def calc_energy(self, label):
        self.set_calculator(label)
        self.energy = self.atoms.get_potential_energy()

    def get_energy(self, label='TiO2'):
        try:
            return self.energy
        except:
            self.calc_energy(label)
            return self.energy

    def constraint_to_string(self):
        a = '{ \n'
        for jj in self.atoms.constraints[0].index:
            a += '{} 1.0 0 0 \n'.format(jj)
            a += '{} 1.0 1 0 \n'.format(jj)
            a += '{} 1.0 0 1 \n'.format(jj)
        a += '}'
        return a











1