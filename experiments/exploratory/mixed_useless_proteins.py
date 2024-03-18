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
frac_useful = [[1]]#, i] for i in np.linspace(0.5, 1, num_species)]
metabolic_rates = [[20.0]]#, 1.0]]# for _ in range(num_species)]
hierarchy = [[0]]#, 1]]# for _ in range(num_species)]
strategy = ['dynamic']# for _ in range(num_species)]
K = [[1E-5]]# 1E-7]]# for _ in range(num_species)]
n = [[2]]#, 1]]# for _ in range(num_species)]
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
nutrients = {'init_concs': [0.1]}#, 0.03]}
ecosystem = diaux.model.Ecosystem(bugs, nutrients)
ecosystem.preculture()#steadystate=False)

#%%
species_df, nutrient_df = ecosystem.grow(10, dt=1E-3, bottleneck={'type':'time', 'interval':2, 'target':0.04})

#%%
plt.plot(species_df['time_hr'], species_df['phi_Rb'], '-', lw=1)
#%%
fig, ax = plt.subplots(1,2, figsize=(6, 4))
grouped = species_df.groupby(['time_hr'])[['M']].sum().reset_index()
ax[0].plot(grouped['time_hr'], grouped['M']/diaux.model.OD_CONV, 'k-', lw=1.5)
# ax[1].set_ylim([1E-2, 1])
# ax[1].set_xlim([0, 25])
# ax[1].set_ylim([1E-5, 1E-2])
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)

ax[0].set_ylabel('approximate OD', fontsize=6)
ax[1].set_ylabel('charged tRNA / Km', fontsize=6)


# for a in ax:
#     # a.set_yscale('log')
for g, d in species_df.groupby('species_label'):
    ax[0].plot(d['time_hr'], d['M'] / diaux.model.OD_CONV, lw=1.5)
    ax[1].plot(d['time_hr'], d['tRNA_u'] / bugs[1], lw=1.5)
    # ax.plot(d['time_hr'], d['M_Mb_2'] / d['M'], '--', lw=1)
    # ax.plot