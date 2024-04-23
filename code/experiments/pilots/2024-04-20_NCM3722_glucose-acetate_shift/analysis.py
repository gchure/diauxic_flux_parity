#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.model 
import diaux.viz 
import scipy.signal
cor, pal = diaux.viz.matplotlib_style()
OD_BOUNDS = [5E-2, 1]
DATE = '2024-04-20'
data = pd.read_csv(f'./processed/{DATE}_complete_growth_curves.csv')

# Restrict the data to the valid bounds
data = data[(data['od650nm'] >= OD_BOUNDS[0]) & 
            (data['od650nm'] <= OD_BOUNDS[1])]

# Set up a figure showing all of the replicates
N_COL = 2
N_REP = data['tech_rep'].max()
N_ROWS = int(np.ceil(N_REP / N_COL))

fig, ax = plt.subplots(N_ROWS, N_COL, figsize=(6, 2 * N_ROWS))
ax = ax.ravel()
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel('log optical density', fontsize=6)
for g, d in data.groupby('tech_rep'):
    ax[g-1].plot(d['time_hr'], np.log(d['od650nm']), '.', ms=4, markeredgewidth=0)
    ax[g-1].set_title(f'Device: M{g-1}')
plt.tight_layout()
plt.savefig('./growth_curves.png', bbox_inches='tight', dpi=150)

#%%
VALID_REPS = [1, 2, 3, 4]
data = data[data['tech_rep'].isin(VALID_REPS)]
data.to_csv(f'./processed/{DATE}_valid_shifts.csv', index=False)