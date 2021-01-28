#!/usr/bin/env python3
# -*- coding: UTF-8 no BOM -*-
import damask
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

r = damask.Result('test_shearZX.hdf5')
r.add_Cauchy()
r.add_strain_tensor()
#r.add_Mises('P')
#r.add_Mises('sigma')
#r.add_Mises('epsilon_V^0.0(F)')
mises_cauchy_stress = []
mises_piola_stress = [] 
mises_strain = []
piola = []
strain = []
increment = []
for inc in r.iterate('increments'):
  piola.append(np.average(r.read_dataset(r.get_dataset_location('P'))[:,0,2]))
  strain.append(np.average(r.read_dataset(r.get_dataset_location('epsilon_V^0.0(F)'))[:,0,2])) 
  increment.append(inc)

#for inc in r.iterate('increments'):
# mises_cauchy_stress_field = r.read_dataset(r.get_dataset_location('sigma_vM'))
# mises_piola_stress_field = r.read_dataset(r.get_dataset_location('P_vM'))
# mises_strain_field = r.read_dataset(r.get_dataset_location('epsilon_V^0.0(F)_vM'))
# mises_cauchy_stress.append(np.average(mises_cauchy_stress_field[:]))
# mises_piola_stress.append(np.average(mises_piola_stress_field[:]))
# mises_strain.append(np.average(mises_strain_field[:]))


plt.xlabel('strain')
#plt.plot(mises_strain,mises_cauchy_stress,'r--', mises_strain,mises_piola_stress, 'b:')
plt.plot(strain,piola)
plt.savefig('stress_strain.png')




