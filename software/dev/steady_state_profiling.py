#%%
import numpy as np 
import pandas as pd
import diaux.viz 
import diaux.model
import diaux.quant
import matplotlib.pyplot as plt
import importlib
importlib.reload(diaux.quant)
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()

suballocation = {'strategy': 'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'nu_max': [10, 0.5]}
nutrients = {'init_concs': [0.0005, 1000.0]}
species = diaux.model.FluxParityAllocator(suballocation)
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(od_init = 0.001, equil_time=100)#, init_conc_override=[1.0, 0.0])
species_df, nutrient_df = ecosystem.grow(50, dt=1/60)
steadystates = diaux.quant.profile_steady_states(species_df)


#%%
# ss = pd.read_csv('../../data/simulations/preshift_postshift_sweep_steady_states.csv')
# traj = pd.read_csv('../../data/simulations/preshift_postshift_sweep_species_trajectories.csv')
# ss.head()
#%%
fig, ax = plt.subplots(2,2)
ax = ax.ravel()
ax[0].set_yscale('log')
ax[0].plot(species_df['time_hr'], species_df['M']/diaux.model.OD_CONV, 'k-', lw=1)

ax[1].plot(species_df['time_hr'], species_df['gamma'] * species_df['ribosome_content'], '-', 
           color=cor['primary_gold'], lw=1)

ax[2].plot(species_df['time_hr'], species_df['alloc_stability'], '-',
           color=cor['primary_red'], lw=1)
# ax[2].set_ylim([-0.01, 0.01])
ax[3].plot(species_df['time_hr'], 1 - species_df['tRNA_c_stability'], '-', color=cor['primary_blue'], lw=1)

# ax[2].hlines(1E-5, 0, 15, color='k', lw=2)
# ax[2].set_ylim([1E-6, 1])
# ax[2].set_yscale('log')
ax[3].set_ylim([-1E-2, 1E-2])
for g, d in steadystates.groupby('total_ss_idx'):
    for a in ax:
        a.fill_betweenx(a.get_ylim(), d['time_start'], d['time_stop'],
                       color=cor['pale_black'], alpha=0.5)

#%%
fig, ax = plt.subplots(1, 1)
dt = species_df['time_hr'].values[:-1]
alloc_stability =  np.diff(species_df['alloc_stability'])
lam = np.diff(species_df['gamma'] * species_df['ribosome_content']) 
# ax.plot(dt, np.abs(alloc_stability/species_df['alloc_stability'].values[:-1]), color=cor['primary_blue'], lw=1)
# ax2 = ax.twinx()
ax.plot(dt, np.abs(lam/(species_df['gamma'].values[:-1] * species_df['ribosome_content'].values[:-1])), color=cor['primary_red'], lw=1)
ax.set_yscale('log')
ax.hlines(5E-5, 0, 25)
ax.set_ylim(1E-8, 1E-2)

# ax.set_hlines(0, 25)
# ax2.set_ylim(0, 1E-3)
