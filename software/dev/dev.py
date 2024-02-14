#%%

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import diaux.model
import diaux.viz
import diaux.callbacks
import importlib
importlib.reload(diaux.model)
importlib.reload(diaux.callbacks)
cor, pal = diaux.viz.matplotlib_style()
suballocation = {'strategy': 'static',
                'alpha': [0.5, 0.5],
                'nu_max': [10, 4]}
bugs1 = diaux.model.FluxParityAllocator(suballocation, label=1)
suballocation = {'strategy': 'dynamic',
                'K': [1E-5, 1E-5],
                'n': [1, 1],
                'nu_max': [10, 4]}
bugs2 = diaux.model.FluxParityAllocator(suballocation, label=2)

community = [bugs2, bugs1]
nutrients = {'init_concs': [0.001, 0.001]}
ecosystem = diaux.model.Ecosystem(community, nutrients)
ecosystem.seed()
species, nuts = ecosystem.grow(time=70, 
                               bottleneck={'type':'time', 
                                            'interval':3, 
                                            'target':0.04},
                               extinction_thresh=0.3)

#%%
for g, d in species.groupby('species_label'):
    plt.plot(d['time_hr'], d['M']/1.5E17)
# plt.yscale('log')


#%%

for g, d in species.groupby('species_label'):
    plt.plot(d['time_hr'], d['frequency'])
