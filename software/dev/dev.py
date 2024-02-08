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

suballocation = {'strategy': 'dynamic',
                'K': [1E-5, 1E-5],
                'n': [1, 1],
                'nu_max': [10, 4]}
bugs = diaux.model.FluxParityAllocator(suballocation, label=1)
community = [bugs]
nutrients = {'init_concs': [0.001, 0.001]}
ecosystem = diaux.model.Ecosystem(community, nutrients)


#%%
ecosystem.initialize()
sol = ecosystem.integrate()
