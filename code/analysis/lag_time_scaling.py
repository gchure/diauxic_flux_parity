#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model
import diaux.viz 
import seaborn as sns
cor, pal = diaux.viz.matplotlib_style()

# Load the datasets and do the pruning as necessary
lag_times = pd.read_csv('../../data/simulations/preshift_postshift_sweep_lag_times.csv')
lag_times.dropna(inplace=True) # Removes an edge case 
traj = pd.read_csv('../../data/simulations/preshift_postshift_sweep_species_trajectories.csv')
traj = traj[traj['nu_max_preshift'].isin([1, 5,  10, 15])]


# Select a representative preshift and assign postshift colors
traj_mid = traj[traj['nu_max_preshift']==10]
colors = sns.color_palette('crest', n_colors=len(traj_mid['nu_max_postshift'].unique()))
cmap = {g:c for g, c in zip(traj_mid['nu_max_postshift'].unique(), colors)}
# traj_mid['color'] = [cmap[g] for g in traj_mid['nu_max_postshift'].values]

# Set up the figure canvas, format, and assign colors
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax = ax.ravel()
ax[0].set_yscale('log')
ax[0].set_xlim([0, 9])
ax[0].set_ylim([0.001,50])
for g, d in traj_mid.groupby('nu_max_postshift'):
    ax[0].plot(d['time_hr'], d['M']/diaux.model.OD_CONV, '-', lw=1,
               color=cmap[g])
    ax[1].plot(d['time_hr'], d['gamma'] * d['ribosome_content'], '-', lw=1,
               color=cmap[g])
    ax[2].plot(d['time_hr'], d['gamma']/d['gamma'].values[0], '-',
               lw=1, color=cmap[g])
    ax[3].plot(d['time_hr'].values[1:], np.diff(d['tRNA_c'])/d['tRNA_c'].values[0], '-',
               lw=1, color=cmap[g])
# ax[3].set_yscale('log')
# ax[3].set_ylim([1E-2,2])
ax[3].set_ylim([-0.05, 0.05])
ax[3].set_xlim(2.5, 5)
# fig, ax = plt.subplots(1,1)
# ax.plot(lag_times['nu_max_preshift']/lag_times['nu_max_postshift'], 1/ lag_times['lag_time_hr'], 'o')

#%%
