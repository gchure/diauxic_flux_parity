#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model
import diaux.viz 
import seaborn as sns
cor, pal = diaux.viz.matplotlib_style()

# Load the trajectories
traj = pd.read_csv('../../data/simulations/preshift_variation_trajectories.csv')
lag = pd.read_csv('../../data/simulations/preshift_variation_lagtimes.csv')
steadystates = pd.read_csv('../../data/simulations/preshift_variation_steadystates.csv')

for g, d in traj.groupby(['nu_max_preshift']):
    plt.plot(d['time_hr'], d['M']/diaux.model.OD_CONV)
    plt.yscale('log')
    # plt.xlim(0, 15)
    # plt.ylim([0, 1])

#%%
 
_lag = lag[lag['nu_max_postshift'].isin([2])]
for g, d in _lag.groupby(['nu_max_postshift']):
    plt.plot(d['preshift_growth_rate_hr'], d['lag_time_hr'], 'o', label=g)
plt.legend()

#%%
plt.plot(lag['postshift_growth_rate_hr'], lag['lag_time_hr'], 'o')