#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model
import diaux.viz 
import seaborn as sns
cor, pal = diaux.viz.matplotlib_style()

# Load the trajectories
traj = pd.read_csv('../../data/simulations/nu_sweep_trajectories.csv')
lag = pd.read_csv('../../data/simulations/nu_sweep_lagtimes.csv')
steadystates = pd.read_csv('../../data/simulations/nu_sweep_steadystates.csv')

# Choose an example shift
nu_pre = 10
growth_curve = traj[traj['nu_max_preshift']==nu_pre]

#%%
# Set up a figure canvas showing the effect of the growth behavior
cmap = sns.color_palette('crest', n_colors=len(growth_curve['nu_max_postshift'].unique()))
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_ylabel('approximate optical density [OD$_{600nm}$ / mL]', fontsize=6)
ax.set_yscale('log')
ax.set_xlim([0, 6])
ax.set_ylim([0.04, 10]) 
for i, (g, d) in enumerate(growth_curve.groupby('nu_max_postshift')):
    _d = d[d['M'] >= lag['shift_biomass'].values[0]]
    __d = d[d['M'] <= lag['shift_biomass'].values[0]]
    ax.plot()
    ax.plot(_d['time_hr'], _d['M']/diaux.model.OD_CONV, lw=1, color=cmap[i])
    ax.plot(__d['time_hr'], __d['M']/diaux.model.OD_CONV, lw=1, color=cor['primary_red'])
plt.savefig('../../figures/sandbox/changing_shift_magnitude.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))

cmap = sns.color_palette('crest', n_colors=len(lag['nu_max_preshift'].unique()))
for i, (g, d) in  enumerate(lag.groupby('nu_max_preshift')):
    ax.plot(d['preshift_growth_rate_hr'] / d['postshift_growth_rate_hr'], d['lag_time_hr'], '-',
            color=cmap[i]) 
    ax.plot(d['preshift_growth_rate_hr'] / d['postshift_growth_rate_hr'], d['lag_time_hr'], 'o',
            markeredgecolor='k', markeredgewidth=0.5, markerfacecolor=cmap[i], ms=4)
#%%
colors = [cor['primary_black'], cor['primary_blue'], cor['primary_red'],
          cor['primary_green'], cor['primary']]
# fig, ax = plt.subplots
#%%
# for i, (g, d) in enumerate(lag.groupby('nu_max_preshift')):
    # plt.plot(d['postshift_growth_rate_hr'],
            #  1/d['lag_time_hr'],
            #  'o') 
for g, d in lag.groupby(['nu_max_postshift']):
    plt.plot(d['preshift_growth_rate_hr'], 1/d['lag_time_hr'], '-o')
    # break
# plt.yscale('log')
# plt.xscale('log')