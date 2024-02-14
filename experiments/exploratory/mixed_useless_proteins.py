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
num_species = 8
frac_useful = [[1, i] for i in np.linspace(0.5, 1, num_species)]
metabolic_rates = [[10.0, 2.0] for _ in range(num_species)]
hierarchy = [[0, 1] for _ in range(num_species)]
strategy = ['dynamic' for _ in range(num_species)]
K = [[1E-7, 1E-7] for _ in range(num_species)]
n = [[2, 2] for _ in range(num_species)]

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

#%%
# With the species in place, set up the ecosystem with defined nutrient parameters
nutrients = {'init_concs': [0.0005, 0.001]}
ecosystem = diaux.model.Ecosystem(bugs, nutrients)
ecosystem.seed()
species_df, nutrient_df = ecosystem.grow(80, bottleneck={'type':'time', 'interval':20, 'target':0.001}, 
                                         dt=None)

#%%
fig, ax = plt.subplots(2,2, figsize=(6, 4))
grouped = species_df.groupby(['time_hr'])[['M']].sum().reset_index()
ax[0,0].plot(grouped['time_hr'], grouped['M']/diaux.model.OD_CONV, 'k-', lw=1.5)
ax[0,0].set_yscale('log')
for g, d in species_df.groupby('species_label'):
    ax[0,0].plot(d['time_hr'], d['M'] / diaux.model.OD_CONV, lw=1.5)
    ax[1, 0].plot(d['time_hr'], d['frequency'])
    # ax.plot(d['time_hr'], d['M_Mb_2'] / d['M'], '--', lw=1)
    # ax.plot(d['time_hr'], d['phi_Mb_2'], 'k-')
    # ax.plot(d['time_hr'], d['M_Rb'] / d['M'], 'b:', lw=1)
    # ax.plot(d['time_hr'], d['tRNA_c'] / d['tRNA_u'])



#%%
