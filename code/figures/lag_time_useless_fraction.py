#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import diaux.model 
import diaux.viz 
import seaborn as sns
cor, pal = diaux.viz.matplotlib_style()

trajectories = pd.read_csv('../../data/simulations/useful_frac_sweep_trajectories.csv')
steadystates = pd.read_csv('../../data/simulations/useful_frac_sweep_steadystates.csv')
lagtimes = pd.read_csv('../../data/simulations/useful_frac_sweep_lagtimes.csv')

n_useless = len(lagtimes['useful_fraction_post'].unique())
cmap = sns.color_palette("flare", n_colors=n_useless)
# fracs = lagtimes['useful_fraction_post'].unique()[::8]

fig, ax = plt.subplots(1, 1, figsize=(3, 2))
for i, (g, d) in enumerate(trajectories[trajectories['nu_max_postshift']==trajectories['nu_max_postshift'].unique()[6]].groupby('useful_fraction_post')):
    _d = d[d['M'] <= lagtimes['shift_biomass'].values[0]]
    __d = d[d['M'] > lagtimes['shift_biomass'].values[0]]
    ax.plot(__d['time_hr'], __d['M']/diaux.model.OD_CONV, '-', color=cmap[i], lw=1)
    ax.plot(_d['time_hr'], _d['M']/diaux.model.OD_CONV, '-', color=cor['primary_black'], lw=1)
ax.set_ylim([0.04, 10])
ax.set_xlim([0, 10])
ax.set_yscale('log')
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_ylabel('approximate optical density [OD$_{600nm}$ / mL]', fontsize=6)
plt.savefig('../../figures/presentations/useless_protein_growth_curves.pdf')


#%%
fig, ax = plt.subplots(1,1, figsize=(2,2))

markers = ['s', 'o', '^', 'v', '>', '<', 'X', 'D', 'h', 'p']
cmap = sns.color_palette('flare', n_colors=10)
for i, (g, d) in enumerate(lagtimes[lagtimes['nu_max_postshift'].isin(lagtimes['nu_max_postshift'].unique()[::3])].groupby('nu_max_postshift')):
    for j, (_g, _d) in enumerate(d.groupby('useful_fraction_post')):
        ax.plot((1 - _d['useful_fraction_post']), 1/_d['lag_time_hr'], 
            markers[i], color=cmap[j], ms=4, markeredgewidth=0.5,
            markeredgecolor=cor['primary_black'], label='__nolegend__')
    ax.plot([], [], marker=markers[i], linestyle='none', color=cor['light_black'],
            markeredgewidth=0.5, markeredgecolor=cor['primary_black'],
            ms=4, label=f'{g:0.1f}')

l = ax.legend(fontsize=5, title=r'$\nu_{max}$ [hr$^{-1}$]' + '\npostshift')
l.get_title().set_fontsize(5)
ax.set_xlabel('useless protein fraction\n$f_x$', fontsize=6)
ax.set_ylabel('1/$\Delta$ [hr$^{-1}$]\ninverse lag time', fontsize=6)
plt.savefig('../../figures/presentations/useless_proteins_lag_times.pdf')