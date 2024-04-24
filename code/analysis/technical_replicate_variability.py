#%%
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import diaux.viz 
import seaborn as sns
cor, pal = diaux.viz.matplotlib_style()


growth_curves = pd.read_csv('../../data/bioreactor/pilot_diauxic_shifts.csv')
growth_curves['idx'] = growth_curves.groupby(['date', 'tech_rep']).ngroup() + 1
cmap = sns.color_palette('husl',n_colors=10)
# Plot all of the curves together
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_xlabel('time from shift [hr]', fontsize=6)
ax.set_ylabel('log relative optical density [a.u.]', fontsize=6)
for g, d in growth_curves.groupby('idx'):
    ax.plot(d['time_hr_centered'], np.log(d['od650nm_centered']), '-', label=g,
            color=cmap[g-1])
leg = ax.legend(title='replicate', fontsize=6)
leg.get_title().set_fontsize(5)
ax.set_title('replicate variability: glucose-acetate shift', fontsize=6)
plt.savefig('../../figures/sandbox/glucose_acetate_replicate_variability.pdf')


lags = pd.read_csv('../../data/bioreactor/pilot_lag_times.csv')
fig, ax = plt.subplots(1,1, figsize=(2,2))
ax.hist(lags['lag_time'])
ax.set_xlabel('lag time [hr]', fontsize=6)
ax.set_ylabel('count', fontsize=6)
#%%
fig, ax = plt.subplots(1,1, figsize=(3,2))
for g, d in growth_curves.groupby('idx'):
    
    ax.plot(d['time_hr_centered'].values[1:],
            np.diff(d['log_od650nm_filtered'])/np.diff(d['time_hr']),
            '-', label=g, color=cmap[g-1])
leg = ax.legend(title='replicate', bbox_to_anchor=(1,1))
leg.get_title().set_fontsize(6)
ax.set_ylim([-0.3, 1.2])
ax.set_ylabel('instantaneous growth rate [hr$^{-1}$]', fontsize=6)
ax.set_xlabel('time from shift [hr]', fontsize=6)
plt.savefig('../../figures/sandbox/glucose_acetate_instantaneous_growth_rate.pdf')