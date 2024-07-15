#%%
import importlib
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model 
import diaux.viz 
import diaux.callbacks
importlib.reload(diaux.callbacks)
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()

suballocation  = {'strategy': 'static',
                  'alpha': [1],
                  'nu_max': [10]}

species = diaux.model.FluxParityAllocator(suballocation)
nutrients = {'init_concs': [1]}
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(equil_time=100, verbose=False)

bottleneck = {'type':'biomass',
              'threshold': 1E18,
              'target': 1E16,
            #   'nutrient_carryover': True
              }
species_df, nutrient_df = ecosystem.grow(20, dt=1/60, bottleneck=bottleneck)

plt.plot(species_df['time_hr'], species_df['M'])
# plt.plot(nutrient_df['time_hr'], nutrient_df['env_conc'])

#%%
# Long term competition
static_suballocation  = {'strategy': 'static',
                  'alpha': [0.6, 0.4],
                  'nu_max': [3, 1],
                  'Y': [3E19, 1E19]}
dynamic_suballocation = {'strategy': 'hierarchical',
                         'K': [1E-5, 1E-5],
                         'n': [1, 1],
                         'nu_max': [3, 1],
                         'Y':[3E19, 1E19]}

static_species = diaux.model.FluxParityAllocator(static_suballocation,
                                                 label=0)
dynamic_species = diaux.model.FluxParityAllocator(dynamic_suballocation,
                                                  label=1)
bottleneck = {'type':'biomass',
              'threshold': 1E18,
              'target': 1E16,
              }
nutrients = {'init_concs': [0.001, 1]}
ecosystem = diaux.model.Ecosystem([static_species, dynamic_species], nutrients)
ecosystem.preculture(equil_time=100)
#%%
species_df, nutrient_df = ecosystem.grow(100, dt=1/60, 
                                         bottleneck=bottleneck)
for g, d in species_df.groupby('species_label'):
    plt.plot(d['time_hr'], d['M'], label=g)

grouped = species_df.groupby('time_hr').sum().reset_index()
plt.plot(grouped['time_hr'], grouped['M'], '--', label='total')
plt.yscale('log')

#%%
n_species = 3 
Ks = 1E-5  * np.ones(n_species)
species = []
suballocation = {'strategy': 'hierarchical',
                 'nu_max': [3, 1],
                 'K': [1E-5, 1E-5],
                 'Y': [3E19, 1E19],
                 'n': [1, 1]}
for i, k in enumerate(Ks):
    # suballocation['K'] = [1E-5, k] 
    _species = diaux.model.FluxParityAllocator(suballocation,
                                               label=i)
    species.append(_species)
nutrients = {'init_concs': [10.0,  0]}
ecosystem = diaux.model.Ecosystem(species, nutrients)
ecosystem.preculture(equil_time=100)

#%%
bottleneck = {'type':'biomass',
              'threshold': 5E16,
              'target': 1E15,
              }
species_df, nutrient_df = ecosystem.grow(1, dt=1/60)
             
#%%
for g, d in species_df.groupby('species_label'):
    plt.plot(d['time_hr'], d['M'], label=g)
    # plt.plot(d['time_hr'], d['frequency'], label=g)
plt.yscale('log')    
plt.legend()
# plt.legend()