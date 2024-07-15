#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import diaux.viz 
import diaux.model
cor, pal = diaux.viz.matplotlib_style()

#%%
# Set the species 
suballocation = {'strategy': 'static',
                 'nu_max': [10.0, 3.0],
                 'alpha': [0.6, 0.4]}
nutrients = {'init_concs': [0.001, 10]}
species = diaux.model.FluxParityAllocator(suballocation, label=1)
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(equil_time=100)
species_df, nut_df = ecosystem.grow(20, dt=1/60)

#%%
plt.plot(species_df['time_hr'], species_df['gamma'] * species_df['ribosome_content'])
plt.ylim([0, 0.5])