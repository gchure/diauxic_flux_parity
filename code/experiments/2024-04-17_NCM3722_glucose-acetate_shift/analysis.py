#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.model 
import diaux.viz 
import scipy.signal
cor, pal = diaux.viz.matplotlib_style()
OD_BOUNDS = [5E-2, 1]
WINDOW = [5]
data = pd.read_csv('./processed/2024-04-17_growth_curves.csv')

# Restrict the data to the valid bounds
data = data[(data['od650nm'] >= OD_BOUNDS[0]) & 
            (data['od650nm'] <= OD_BOUNDS[1])]

# Play with the first replicate
rep = data[data['replicate']==1]

plt.plot(rep['time_hr'], rep['od650nm'], '-', lw=1)
plt.yscale('log')
#%%

fig, ax = plt.subplots(1,1, figsize=(4,2.5))
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_ylabel('OD$_{650nm}$', fontsize=6)
for g, d in data.groupby('replicate'):
    ax.plot(d['time_hr'].values, d['od650nm'], label=g)
ax.set_yscale('log')
ax.set_ylim([5E-2, 2])

ax.set_xlim([1, 8])
ax.legend(title='replicate')

# ax.set_ylim([0, 0.1])
#%%
scipy.signal.medfilt(d['od650nm'].values, 7)
