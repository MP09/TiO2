import numpy as np
import sys, os
from basin_hoping import basin_hopping, TiO2
from pertubations import pertubation_dict
from timeit import default_timer as dt

# Arrange inputs:
ncalc = int(sys.argv[1])                 # How many calculations to do.
ninit = int(sys.argv[2])                 # What number is the first calculation?
pertubation_number = int(sys.argv[3])    # Which pertubation to use.
ml_setting = int(sys.argv[4])            # Whether or not to use autobag.
max_iterations = int(sys.argv[5])        # Maximum number of iterations before terminating BH

# Pertubation settings:
pertubation_info = pertubation_dict[pertubation_number]
pertubation_full_name = pertubation_info[0]
pertubation_short_name = pertubation_info[1]
pertubation = pertubation_info[2]

# Autobag settings: 
if ml_setting == 0:
    use_ML = True
    ML_name = 'ML_'
else:
    use_ML = False
    ML_name = ''

# Naming:
traj_name = ML_name+pertubation_short_name+'{}.traj'
npy_name = ML_name+pertubation_short_name+'{}.npy'

# Start up for loop:
for ii in range(ninit, ninit+ncalc):
    try: 
        # Intialize basin hopping:
        Ba = basin_hopping(traj_name=traj_name.format(ii), use_ML=use_ML)
        Ba.max_iterations = max_iterations
        Ba.set_perturbation(pertubation)

        t0 = dt()
        conv, niter, E = Ba.start()
        time = dt()-t0
        result = np.array([conv, niter, E, time])
        np.save(npy_name.format(ii), result)

    except FileNotFoundError:
        os.rename('dftb.out', ML_name+'dftb_{}.fail'.format(ii))





    
