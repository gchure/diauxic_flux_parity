#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import diaux.model
import diaux.viz
import scipy.integrate
import importlib
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()

suballocation = {'strategy': 'static',
                'alpha': [0.8, 0.2],
                'nu_max': [10, 4]}
bugs1 = diaux.model.FluxParityAllocator(suballocation, label=1)
suballocation = {'strategy': 'dynamic',
                'K': [1E-5, 1E-5],
                'n': [1, 1],
                'nu_max': [10, 4]}
bugs2 = diaux.model.FluxParityAllocator(suballocation, label=2)

community = [bugs1, bugs2]
nutrients = {'init_concs': [0.001, 0.001]}
ecosystem = diaux.model.Ecosystem(community, nutrients)

#%%
ecosystem.seed()
species, nuts = ecosystem.grow(time=10)
# species, nut = ecosystem._integrate(time_range=[0, 10], 
                        #   solver_kwargs={'t_eval':np.linspace(0, 10, 200)})

#%%
for g, d in species.groupby('species_label'):
    plt.plot(d['time_hr'], d['frequency'])
