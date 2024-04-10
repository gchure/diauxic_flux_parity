#%%
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
import pandas as pd 
import diaux.model
import diaux.viz 
import seaborn as sns
cor, pal = diaux.viz.matplotlib_style()


# %% Set up two examples
nu_max = [10.0, 5.0]
static_suballocation = {'strategy':'static',
                        'alpha':[0.7, 0.3],
                        'K': [1E-5, 1E-5],
                        'n': [1, 1],
                        'nu_max':nu_max}

dynamic_suballocation = {'strategy':'hierarchical',
                         'nu_max':nu_max,
                         'K': [1E-5, 1E-5],
                         'n': [1, 1]}

prop_suballocation = {'strategy':'proportional',
                      'nu_max':nu_max,
                       'K': [1E-5, 1E-5],
                       'n': [1, 1]}

static_FPA = diaux.model.FluxParityAllocator(static_suballocation)

dynamic_FPA = diaux.model.FluxParityAllocator(dynamic_suballocation)

prop_FPA = diaux.model.FluxParityAllocator(prop_suballocation)


#%%
# Define the ecosystems
nutrients = {'init_concs':[0.001, 0.001]}
static_eco = diaux.model.Ecosystem([static_FPA], nutrients)
dynamic_eco = diaux.model.Ecosystem([dynamic_FPA], nutrients)
prop_eco = diaux.model.Ecosystem([prop_FPA], nutrients)

#%% Preculture the two systems to their initial steadystates
static_eco.preculture(od_init=0.05)
dynamic_eco.preculture(od_init=0.05)
prop_eco.preculture(od_init=0.05)

#%% Growt the communities
static_species, static_nuts = static_eco.grow(3, dt=0.001)
static_species['inst_lam'] = static_species['gamma'] * static_species['ribosome_content']
dynamic_species, dynamic_nuts = dynamic_eco.grow(3, dt=0.001)
dynamic_species['inst_lam'] = dynamic_species['gamma'] * dynamic_species['ribosome_content']
prop_species, prop_nuts = prop_eco.grow(3, dt=0.001)
prop_species['inst_lam'] = prop_species['gamma'] * prop_species['ribosome_content']

#%% Set up the axes
gs = GridSpec(4, 3)
fig = plt.figure(figsize=(6, 5))
biomass_ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]
growth_rate_ax = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])]
total_allocation_ax = [fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1]), fig.add_subplot(gs[2, 2])] 
metabolic_allocation_ax = [fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1]), fig.add_subplot(gs[3, 2])] 

for ta, ma in zip(total_allocation_ax, metabolic_allocation_ax):
    for a in [ta, ma]:
        a.set_facecolor('none')  
        a.grid(False)
        a.set_yticklabels([])
        a.set_xticks([0, 0.5,  1])
        a.set_xticklabels(['0', '0.5', '1'])
        a.set_ylim([0.5, 4])
        a.set_xlim([-0.1, 1])
    for i, l in enumerate(['$t_1$', '$t_2$', '$t_3$']):
        for a in [ta,ma]:
            a.plot(-0.05, 3-i, 'o', ms=8, markeredgecolor=cor['primary_black'],
                markerfacecolor='w', markeredgewidth=0.75)
            a.text(-0.078, 3-i-0.11, l, fontsize=6, fontweight='bold', color=cor['primary_black'])
for a in growth_rate_ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel(r'$\lambda(t)/\lambda_0$' + '\n relative growth rate', fontsize=6)
    a.set_yticks([0, 0.5, 1.0])
    a.set_ylim([-0.05, 1.1])
for a in biomass_ax:
    a.set_ylim([0, 10])
    a.set_ylabel('$M(t)/M_0$\nrelative biomass', fontsize=6)
    a.set_xticklabels([])

# Plot the time points to demonstrate the allocation behavior
static_timepoints = [0.75, 1.165, 1.9]
dynamic_timepoints = [0.5, 0.99, 2.3]
prop_timepoints = [0.5, 1.4, 1.655]

for i, st, pt, dt, l in zip([3, 2, 1], static_timepoints, prop_timepoints, dynamic_timepoints, ['$t_1$', '$t_2$', '$t_3$']):
    _static = static_species[static_species['time_hr']>=st].iloc[0]
    _dynamic = dynamic_species[dynamic_species['time_hr']>=dt].iloc[0]
    _prop = prop_species[prop_species['time_hr']>=pt].iloc[0]
    for j, _d, d in zip([0, 1, 2], [_static, _prop, _dynamic], [static_species, prop_species, dynamic_species]): 
        # Plot the points to indicate the time points and label them
        biomass_ax[j].plot(_d['time_hr'], _d['M']/d['M'].values[0],
                'o', ms=8, markeredgewidth=0.75, markerfacecolor='w',
                markeredgecolor=cor['primary_black'],zorder=1000) 
        biomass_ax[j].text(_d['time_hr']-0.07, _d['M']/d['M'].values[0]-0.3,
               l, color=cor['primary_black'], fontsize=6, fontweight='bold',
               zorder=1001)
        growth_rate_ax[j].plot(_d['time_hr'], _d['inst_lam']/d['inst_lam'].values[0],
                'o', ms=8, markeredgewidth=0.75, markerfacecolor='w',
                markeredgecolor=cor['primary_black'],zorder=1000) 
        growth_rate_ax[j].text(_d['time_hr']-0.07, _d['inst_lam']/d['inst_lam'].values[0]-0.03,
               l, color=cor['primary_black'], fontsize=6, fontweight='bold',
               zorder=1001)


        # Plot the total resource allocation
        total_allocation_ax[j].barh(i, 0.55, color=cor['light_black'],
                                    edgecolor='w',
                                    linewidth=0.5)
        total_allocation_ax[j].barh(i, _d['phi_Rb'], left=0.55, color=cor['primary_gold'],
                                    edgecolor='w',linewidth=0.5)
        total_allocation_ax[j].barh(i, 1 - 0.55 - _d['phi_Rb'], 
                                    left=0.55 + _d['phi_Rb'], color=cor['purple'],
                                    edgecolor='w',linewidth=0.5)

        # Plot the metabolic suballocation        
        metabolic_allocation_ax[j].barh(i, _d['alpha_1'], color=cor['primary_red'],
                                        edgecolor='w', linewidth=0.5) 
        metabolic_allocation_ax[j].barh(i, _d['alpha_2'], left=_d['alpha_1'], color=cor['pale_red'],
                                        edgecolor='w', linewidth=0.5) 


# Plot where the nutrients exhaust
for i, n in enumerate([static_nuts, prop_nuts, dynamic_nuts]):
    _n = n[n['env_conc']==0]
    for j, a in enumerate([biomass_ax[i], growth_rate_ax[i]]):
        if j == 0:
            ylims = [1, 15]
        else:
            ylims = [-0.1, 1.1]
        ls = ['-', '--']
        for g, d in _n.groupby('nutrient_label'):
            a.vlines(d['time_hr'].values[0], ylims[0], ylims[1], color=cor['light_black'],
                     lw=1, linestyle=ls[g-1], zorder=1)

        a.set_ylim(*ylims)
# Plot the growth curves
for i, d in enumerate([static_species, prop_species, dynamic_species]):
    d = d.copy()
    biomass_ax[i].plot(d['time_hr'], d['M']/d['M'].values[0],'-', lw=1, 
           color=cor['primary_black'])
    growth_rate_ax[i].plot(d['time_hr'], d['inst_lam'] / d['inst_lam'].values[0],
                           '-', lw=1, color=cor['primary_black'])

plt.savefig('../../figures/doc/static_vs_dynamic_vs_prop_plots.pdf')