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
                'n': [4, 1],
                'nu_max': [10, 4]}
bugs1 = diaux.model.FluxParityAllocator(suballocation, label=1)
suballocation = {'strategy': 'dynamic',
                'K': [5E-6, 1E-5],
                'n': [1, 1],
                'nu_max': [10, 4],
                'hierarchy': [1, 0]}
bugs2 = diaux.model.FluxParityAllocator(suballocation, label=2)

community = [bugs1] #, bugs2]
nutrients = {'init_concs': [0.001, 0.001]}
ecosystem = diaux.model.Ecosystem(community, nutrients)


#%%
ecosystem.initialize()
sol = ecosystem.integrate()#solver_kwargs={'method':'LSODA'})

#%%
y = sol.y.T
M = [np.sum(y[i][:4]) for i in range(len(y))]
plt.plot(sol.t, np.array(M)/diaux.model.OD_CONV)
