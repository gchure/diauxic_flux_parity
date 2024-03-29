#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.viz
import diaux.model
import importlib
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()


# Set up species of suballocators with different fractions of the 
# secondary metabolic class deemed to be useful
num_species = 1
frac_useful = [[1, 1]]#, i] for i in np.linspace(0.5, 1, num_species)]
metabolic_rates = [[20.0, 1.0]]#, 1.0]]# for _ in range(num_species)]
hierarchy = [[0, 1]]#, 1]]# for _ in range(num_species)]
strategy = ['dynamic']# for _ in range(num_species)]
K = [[1E-6, 1E-6]]# 1E-7]]# for _ in range(num_species)]
n = [[2, 2]]#, 1]]# for _ in range(num_species)]
# c = {'Km_c': 1E-3}

bugs = []
for i in range(num_species):
    suballocation = {'strategy':strategy[i],
                     'nu_max': metabolic_rates[i],
                     'frac_useful': frac_useful[i],
                     'hierarchy': hierarchy[i],
                     'K': K[i],
                     'n': n[i]}
    FPA = diaux.model.FluxParityAllocator(suballocation, label=(i+1))
                   
    bugs.append(FPA)
# With the species in place, set up the ecosystem with defined nutrient parameters
nutrients = {'init_concs': [0.001, 0.001]}#, 0.03]}
ecosystem = diaux.model.Ecosystem(bugs, nutrients)
ecosystem.preculture(init_conc_override=[0.1, 0])

#%%
species_df, nutrient_df = ecosystem.grow(50, bottleneck={'type':'time', 'interval':15, 'target':0.04})
#%%
species_df
#%%
plt.plot(species_df['time_hr'], species_df['M']/ diaux.model.OD_CONV, '-', lw=1)
plt.yscale('log')