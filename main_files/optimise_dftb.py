#!/usr/bin/env python
#-*- coding: utf-8 -*-
from __future__ import division

import os
from ase import Atoms
from ase.optimize import QuasiNewton, BFGS 
from ase.ga.relax_attaches import VariansBreak
from ase.calculators.dftb import Dftb
from ase.io import write,read
from ase.visualize import view
from gpaw import *
from gpaw import extra_parameters
from ase.constraints import FixAtoms

inputprefix='nonoptimised'
outputprefix='optimisation_dftb'

t=read(inputprefix+'.traj')

t.set_calculator(Dftb(label='C',
                      atoms=t,
                      Hamiltonian_SCC='No',
                      Hamiltonian_Eigensolver='Standard {}',
                      Hamiltonian_MaxAngularMomentum_='',
                      Hamiltonian_MaxAngularMomentum_C='"p"',
                      Hamiltonian_Charge='0.000000',
                      Hamiltonian_Filling ='Fermi {',
                      Hamiltonian_Filling_empty= 'Temperature [Kelvin] = 0.000000',
                      ))

niter = 0
forcemax = 0.05
while (t.get_forces()**2).sum(axis = 1).max()**0.5 > forcemax and niter < 10:
    try:
        dyn = BFGS(t, trajectory=outputprefix+'.traj', 
                   logfile=outputprefix+'.log')
        vb = VariansBreak(t, dyn)
        dyn.attach(vb)
        dyn.run(fmax=forcemax, steps=1000)
    except Exception, e:
        print e
    niter += 1

#E = t.get_potential_energy()
#F = t.get_forces()

#write(outputprefix+'.traj',t)
