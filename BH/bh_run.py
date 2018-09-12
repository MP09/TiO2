import numpy as np
from timeit import default_timer as dt
from basin_hoping import basin_hopping, TiO2
import sys, os

if sys.argv[3] == '1':
    use_ML = True
else:
    use_ML = False


label = sys.argv[4]
traj_name = 'BH_{}.traj'.format(sys.argv[1])
try: 
    Ba = basin_hopping(traj_name=traj_name, use_ML=use_ML)
    Ba.max_iterations = int(sys.argv[2])
    Ba.set_perturbation(label)

    t0 = dt()
    conv, niter, E = Ba.start()
    time = dt()-t0

    if conv == True:
        a = np.array([1, niter, E, time])
    else:
        a = np.array([0, niter, E, time])

    if use_ML:
        np.save('ML_BH_results_{}.npy'.format(sys.argv[1]), a)
        print('Saved some stuff')
    else:
        np.save('BH_results_{}.npy'.format(sys.argv[1]), a)
        print('Saved some stuff')
except FileNotFoundError:
    if use_ML:
        os.rename('dftb.out', 'ML_dftb_{}.fail'.format(sys.argv[1]))
        os.rename('detailed.out', 'ML_detailed_{}.fail'.format(sys.argv[1]))
        #os.rename('geo_end.xyz', 'ML_geo_end_{}.xyz'.format(sys.argv[1]))
        #os.remove('band.out')
    else:
        os.rename('dftb.out', 'dftb_{}.fail'.format(sys.argv[1]))
        os.rename('detailed.out', 'detailed_{}.fail'.format(sys.argv[1]))
        #os.rename('geo_end.xyz', 'geo_end_{}.xyz'.format(sys.argv[1]))
        #os.remove('band.out')
    





