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

#%%
_lag = lag[lag['nu_max_postshift'].isin([1, 2, 3, 4])]
cmap = sns.color_palette('crest', n_colors=len(lag['nu_max_preshift'].unique()))
markers = ['s', 'o', '^', 'v']
fig, ax = plt.subplots(1, 2, figsize=(4,2))
for i, (g, d) in enumerate(_lag.groupby(['nu_max_postshift'])):
    d = d.copy()
    d['rel_lag'] = d['lag_time_hr'] / d['lag_time_hr'].max()

    for j, a in enumerate(ax):
        if j == 0:
            label = g[0]
        else:
            label = f"{d['postshift_growth_rate_hr'].mean():0.2f}"
        a.plot([], [], linestyle='none', marker=markers[i], color=cor['pale_black'], 
               markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
               label=label, ms=4)
    for j, (_g, _d) in enumerate(d.groupby('nu_max_preshift')):
        ax[0].plot(_d['nu_max_preshift'], _d['rel_lag'], markers[i], 
            ms=4, color=cmap[j], markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
            label='__nolegend__')
        ax[1].plot(_d['preshift_growth_rate_hr'], _d['rel_lag'], markers[i], 
            ms=4, color=cmap[j], markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
            label='__nolegend__')

for i, a in enumerate(ax):
    if i == 0:
        title = r'$\nu_{max}$ [hr$^{-1}$]' +'\npostshift'
    else:
        title = '$\lambda$ [hr$^{-1}$] \npostshift'
    a.set_ylabel('$\Delta/\Delta_{max}$ \n relative lag time', fontsize=6)
    l = a.legend(title=title, handletextpad=0.1, fontsize=5)
    l.get_title().set_fontsize(5)
ax[0].set_xlabel('preshift metabolic rate\n'+r'$\nu_{max}$ [hr$^{-1}$]', fontsize=6)
ax[1].set_xlabel('preshift growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)

plt.tight_layout()
plt.savefig('../../figures/preshift_effects.pdf', bbox_inches='tight')

#%%
_lag = lag[lag['nu_max_preshift'].isin([5, 10, 15, 20])]
cmap = sns.color_palette('crest', n_colors=len(lag['nu_max_postshift'].unique()))
markers = ['s', 'o', '^', 'v']
fig, ax = plt.subplots(1, 2, figsize=(4,2))
for i, (g, d) in enumerate(_lag.groupby(['nu_max_preshift'])):
    d = d.copy()
    d['rel_lag'] = d['lag_time_hr'] / d['lag_time_hr'].max()

    for j, a in enumerate(ax):
        if j == 0:
            label = g[0]
        else:
            label = f"{d['preshift_growth_rate_hr'].mean():0.2f}"
        a.plot([], [], linestyle='none', marker=markers[i], color=cor['pale_black'], 
               markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
               label=label, ms=4)
    for j, (_g, _d) in enumerate(d.groupby('nu_max_postshift')):
        ax[0].plot(_d['nu_max_postshift'], 1/_d['lag_time_hr'], markers[i], 
            ms=4, color=cmap[j], markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
            label='__nolegend__')
        ax[1].plot(_d['postshift_growth_rate_hr'], 1/_d['lag_time_hr'], markers[i], 
            ms=4, color=cmap[j], markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
            label='__nolegend__')

for i, a in enumerate(ax):
    if i == 0:
        title = r'$\nu_{max}$ [hr$^{-1}$]' +'\npreshift'
    else:
        title = '$\lambda$ [hr$^{-1}$] \npreshift'
    a.set_ylabel('$1/\Delta$ [hr$^{-1}$] \n inverse lag time', fontsize=6)
    l = a.legend(title=title, handletextpad=0.1, fontsize=5)
    l.get_title().set_fontsize(5)
ax[0].set_xlabel('postshift metabolic rate\n'+r'$\nu_{max}$ [hr$^{-1}$]', fontsize=6)
ax[1].set_xlabel('postshift growth rate\n$\lambda$ [hr$^{-1}$]', fontsize=6)
plt.tight_layout()
plt.savefig('../../figures/posthift_effects.pdf', bbox_inches='tight')


#%%
fig, ax = plt.subplots(1,1, figsize=(3,2))
cmap = sns.color_palette('crest', n_colors=5)
_traj = traj[(traj['nu_max_postshift']==5) &
             (traj['nu_max_preshift'].isin([5.5, 8, 10, 12, 18]))]
_lag = lag[(lag['nu_max_postshift']==5) & 
           (lag['nu_max_preshift'].isin([5.5, 8, 10, 12, 18]))]
for i, (g, d) in enumerate(_traj.groupby('nu_max_preshift', sort=False)):
    t_hit = _lag[_lag['nu_max_preshift']==g]['shift_time'].values[0]
    ax.plot(d['time_hr'] - t_hit, d['M']/diaux.model.OD_CONV, '-',
            lw=1, color=cmap[i])

ax.set_xlabel('time from shift [hr]', fontsize=6)
ax.set_ylabel('approximate optical density [OD$_{600nm}$ / mL]', fontsize=6)
ax.set_yscale('log')
ax.set_xlim([-1.5, 4])
ax.set_ylim([0.04, 10])