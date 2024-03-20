#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.viz
import diaux.model
import importlib
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()

# Define the species which can eat two nutrients
suballocation = {'strategy': 'dynamic',
                 'nu_max': [10.0, 2.0],
                 'frac_useful': [1, 1],     
                 'K': [1E-5, 1E-5],
                 'n': [1, 1]}
species = diaux.model.FluxParityAllocator(suballocation, label=1)                   

# With the species in place, set up the ecosystem with defined nutrient parameters
# and preculture in the pre-shift condition
nutrients = {'init_concs': [0.001, 0.01]}#, 0.03]}
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(init_conc_override=[0.1, 0])

#%%
# Grow the ecosystem    
species_df, nutrient_df = ecosystem.grow(20, dt=1E-4)

#%% 
# Find the steady-states
alloc_ss = np.abs(1 - species_df['alloc_stability']) <= 0.03
tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= 0.03
tRNA_u_ss = np.abs(1 - species_df['tRNA_u_stability']) <= 0.03
species_df['steady_state'] = alloc_ss #* tRNA_c_ss * tRNA_u_ss
plt.plot(species_df['time_hr'], species_df['steady_state'])
#%%
M_INIT = 0.04 * diaux.model.OD_CONV
M_STAR= (nutrients['init_concs'][0] * species.Y[0] + M_INIT)
lam_1 = species_df.iloc[50]['gamma'] * species_df.iloc[50]['phi_Rb']
t_hit = np.log(M_STAR/M_INIT)/lam_1
t_hit
#%%
plt.hlines(M_STAR/diaux.model.OD_CONV, 0, 20, color=cor['primary_red'],
           linewidth=2)
plt.vlines(t_hit, 0, 1, color=cor['primary_blue'], lw=2)
plt.plot(species_df['time_hr'], species_df['M'] / diaux.model.OD_CONV, lw=2)

# plt.yscale('log')
