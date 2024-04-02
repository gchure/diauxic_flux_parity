#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
import seaborn as sns 
import diaux.model 
import diaux.viz 
cor, pal = diaux.viz.matplotlib_style()

# Define some species
alphas = [[1/3, 2/3], 
          [2/3, 1/3],
          [0.4, 0.6], 
        #   [0.7, 0.3],
          [0.9, 0.1]]
nu_max = [10.0, 10.0]


dil_rate = 0.1 
species = []
for i, a in enumerate(alphas):
    suballocation = {'strategy': 'static',
                     'nu_max':nu_max,
                     'alpha': a,
                     'death_rate':dil_rate}
    FPA = diaux.model.FluxParityAllocator(suballocation,
                                          metabolic_hierarchy=False,
                                          label=i+1)
    species.append(FPA)

#%%

# Set up a static ecosystem
feed_concs = [0.0095, 0.0005]
OD_INIT = 5.89923E16/diaux.model.OD_CONV
inflow_rates = [dil_rate, dil_rate]
degradation_rates = [dil_rate, dil_rate]
nutrients = {'init_concs': feed_concs,
             'feed_concs': feed_concs,
             'inflow_rates': inflow_rates,
             'degradation_rates': degradation_rates}
ecosystem = diaux.model.Ecosystem(species, nutrients)

#%%
ecosystem.preculture(od_init=OD_INIT)
species_df, nutrient_df = ecosystem.grow(200, dt=0.001)

#%%
# Assign colors to the species
# cmap = sns.color_palette('pastel', n_colors=len(species)+2)
colors = [cor['primary_black'], cor['primary_blue'], cor['primary_green'],
          cor['primary_purple'], cor['primary_gold']]
species_df['color'] = [colors[j-1] for j in species_df['species_label'].values]


def instantiate_two_nutrient_figure(comp_ax='freq'):
    """
    Instantiates the default ecosystem figure canvas for growth on two figures. 

    Parameters
    ----------
    comp_ax : str
        The nature of the composition axis. Default option is `freq` which 
        states that the mass frequency of each species in the ecosystem should 
        be shown. If `biomass` is provided, the biomass (in approximate OD units)
        will be generated instead.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure canvas.
    line_ax : matplotlib.axes._axes.Axes
        The axis object which defines the 1-simplex representation of the 
        metabolic strategies and nutrient concentrations.
    bar_ax : matplotlib.axes._axes.Axes
        The axis object which defines the metabolic suballocation representation 
        of the metabolic strategies and nutrient concentrations.
    comp_ax : matplotlib.axes._axes.Axes
        The axis object which defines the ecosystem composition representation
        of the simulation.    
    nut_ax : matplotlib.axes._axes.Axes 
        The axis object which defines the nutrient composition in the ecosystem. 

    """
    fig = plt.figure(figsize=(6,3))#, layout="")
    gs = GridSpec(4, 5)

    # Define the gridspec axes
    line_ax = fig.add_subplot(gs[0, :2])
    bar_ax = fig.add_subplot(gs[1:, :2])
    comp_ax = fig.add_subplot(gs[:2, 2:])
    nut_ax = fig.add_subplot(gs[2:, 2:])

    # Add figure panel titles
    line_ax.set_title('metabolic strategy 1-simplex', fontsize=6, y=1.1)
    bar_ax.set_title('metabolic suballocation', fontsize=6, zorder=1000, y=1.08)
    comp_ax.set_title('ecosystem species composition', fontsize=6)
    nut_ax.set_title('environment nutrient composition', fontsize=6)

    # Add axis labels 
    bar_ax.set_xlabel('species label', fontsize=6)
    bar_ax.set_title('metabolic suballocation', fontsize=6, zorder=1000, y=1.08)
    bar_ax.set_ylabel(r'$\alpha$', fontsize=6)
    if comp_ax == 'freq':
        comp_ax.set_ylabel('mass frequency', fontsize=6)
    elif comp_ax == 'biomass':  
        comp_ax.set_ylabel('approximate OD$_{600nm}$', fontsize=6)
    nut_ax.set_ylabel('concentration [mM]', fontsize=6)
    nut_ax.set_xlabel('time [hr]', fontsize=6)

    # Add default structures to the 1-simplex 
    line_ax.set_facecolor('none')
    line_ax.set_xticks([])
    line_ax.set_yticks([])
    line_ax.set_ylim([-0.1, 0.1])
    return [fig, line_ax, bar_ax, comp_ax, nut_ax]


def set_ecosystem_colors(species_df, 
                         nutrient_df, 
                         species_palette='default',
                         nut_palette='default'):
    """
    Assigns colors to each species and nutrient.

    Parameters
    ----------
    species_df : pandas DataFrame
        A dataframe of the species details as returned by `Ecosystem.grow`.
    nutrient_df : pandas DataFrame
        A dataframe of the nutrient details as returned by `Ecosystem.grow`.
    species_palette : Seaborn color palette or st
        The color palette to be used for the species. Providing 'default' will cycle 
        through the diaux color palette.
    nturient_palette : Seaborn color palette or st
        The color palette to be used for the species. Providing 'default' will cycle 
        through the red hues of the diaux color palette. 

    Returns
    -------
    species_df : pandas DataFrame
        The species dataframe with the added colors.
    nutrient_df : pandas DataFrame
        The nutrient datframe with the added colors.
    """


#%%

    # line_ax.plot([0], [0], 's', markerfacecolor=cor['dark_red'], markeredgecolor=cor['primary_black'])
    # line_ax.hlines(0, 1, '-', lw=0.5, color=cor['light_black'])
    # line_ax.text(-0.0125, -0.012, 'A', fontsize=4, color=cor['pale_red'], fontweight='bold')
    # line_ax.plot([1], [0], 's', markerfacecolor=cor['pale_red'], markeredgecolor=cor['primary_black'])
    # line_ax.text(0.9875,  -0.012, 'B', fontsize=4, color=cor['dark_red'], fontweight='bold')


# Set up a figure canvas
fig = plt.figure(figsize=(6,3))#, layout="")
gs = GridSpec(4, 5)

# Define the gridspec axes
line_ax = fig.add_subplot(gs[0, :2])
bar_ax = fig.add_subplot(gs[1:, :2])
freq_ax = fig.add_subplot(gs[:2, 2:])
nut_ax = fig.add_subplot(gs[2:, 2:])
nut_ax.set_xlabel('time [hr]', fontsize=6)

# Adjust various limits and labels
freq_ax.set_title('ecosystem species composition', fontsize=6)
nut_ax.set_title('environment nutrient composition', fontsize=6)
freq_ax.set_ylabel('mass frequency', fontsize=6)
nut_ax.set_ylabel('concentration [mM]', fontsize=6)
for a in [freq_ax, nut_ax]:
    a.set_yscale('log')
freq_ax.set_xlim(species_df['time_hr'].min(), species_df['time_hr'].max())
freq_ax.set_ylim([1E-3, 2])
freq_ax.set_xticklabels([])
nut_ax.set_xlim(freq_ax.get_xlim()[0], freq_ax.get_xlim()[1])

# Apply necessary boilerplate formatting
line_ax.set_facecolor('none')
line_ax.set_xticks([])
line_ax.set_yticks([])
line_ax.set_ylim([-0.1, 0.1])
line_ax.hlines(0, 1, '-', lw=0.5, color=cor['light_black'])
line_ax.plot([0], [0], 's', markerfacecolor=cor['dark_red'], markeredgecolor=cor['primary_black'])
line_ax.text(-0.0125, -0.012, 'A', fontsize=4, color=cor['pale_red'], fontweight='bold')
line_ax.plot([1], [0], 's', markerfacecolor=cor['pale_red'], markeredgecolor=cor['primary_black'])
line_ax.text(0.9875,  -0.012, 'B', fontsize=4, color=cor['dark_red'], fontweight='bold')
line_ax.set_title('metabolic strategy 1-simplex', fontsize=6, y=1.1)

bar_ax.set_facecolor('w')
# bar_ax.set_xticks(np.arange(0, len(species), 1))
# bar_ax.set_xticklabels([s.label for s in ecosystem.species])
bar_ax.set_xlabel('species label', fontsize=6)
bar_ax.set_title('metabolic suballocation', fontsize=6, zorder=1000, y=1.08)
bar_ax.set_ylabel(r'$\alpha$', fontsize=6)


# Plot the species-level info for the specified index
labels = []
for g, d in species_df.groupby('species_label'):
    labels.append(g)
    tp = d.iloc[-1]
    if tp['frequency'] <= 1E-3:
        marker = 'X' 
    else:
        marker = 'o'
    line_ax.plot(1-tp['alpha_1'], 0, marker, label='__nolegend__', ms=4, color=tp['color'],
                 markeredgecolor=cor['primary_black'],
                 markeredgewidth=0.3)
    bottom = 0
    colors = ['purple', 'light_purple']
    for j in range(2):
        suballoc = 1-d.iloc[-1][f'alpha_{j+1}']
        bar_ax.bar(g, suballoc, bottom=bottom, color=cor[colors[j]],
                   label='__nolegend__')
        bottom += suballoc

    bar_ax.plot(g, -.05, marker=marker, ms=4, color=tp['color'], markeredgecolor=cor['primary_black'],
                markeredgewidth=0.3, label='__nolegend__')
bar_ax.set_xticks(labels)

# Plot the species frequency and nutrient concentrations
for g, d in species_df.groupby(['species_label', 'color']):
    freq_ax.plot(d['time_hr'], d['frequency'], lw=1.5, color=g[1])

nut_color = [cor['dark_red'], cor['primary_red']]
labels = ['A', 'B']
for g, d in nutrient_df.groupby('nutrient_label'):
    nut_ax.plot(d['time_hr'], d['env_conc']*1E3, '-', lw=1.5, color=nut_color[g-1],
                label=labels[g-1])
nut_ax.legend(fontsize=5)
nut_tp = nutrient_df[nutrient_df['time_hr']==nutrient_df['time_hr'].max()].copy()
nut_tp['feed_frac'] = nut_tp['feed_conc'] / np.sum(nut_tp['feed_conc'])
nut_tp['env_frac'] = nut_tp['env_conc'] / np.sum(nut_tp['env_conc'])
_nut_tp = nut_tp[nut_tp['nutrient_label']==1]
line_ax.plot(1-_nut_tp['feed_frac'].values[0], 0.02, 'v', ms=4, 
             markerfacecolor=cor['pale_red'],
             markeredgecolor=cor['primary_black'],
             markeredgewidth=0.5, label='nutrient (feed)')
line_ax.plot(1-_nut_tp['env_frac'].values[0], -0.02, '^', ms=4, 
             markerfacecolor=cor['red'],
             markeredgecolor=cor['primary_black'],
             markeredgewidth=0.5, label='nutrient (env)')
line_ax.plot([], [], 'o',markerfacecolor='w', markeredgecolor='k',
             markeredgewidth=0.3, label='species', ms=4)
line_ax.legend(frameon=False, fontsize=4, ncol=3, bbox_to_anchor=(1.05,1.1))


# Plot the nutrient breakdown
bar_lims = bar_ax.get_xlim()

bar_ax.set_xlim(bar_lims)

bar_ax.hlines(1-_nut_tp['env_frac'].values[0], bar_lims[0], bar_lims[1], lw=1, 
            color=cor['red'], label='nutrient (env)')
bar_ax.hlines(1-_nut_tp['feed_frac'].values[0], bar_lims[0], bar_lims[1], lw=1, 
            color=cor['pale_red'], label='nutrient (feed)')

bar_ax.plot([], [], 's', markerfacecolor=cor['purple'], label=r'$\alpha_A$')
bar_ax.plot([], [], 's', markerfacecolor=cor['pale_purple'], label=r'$\alpha_B$')
bar_ax.legend(frameon=False, ncol=4, handlelength=0.5, fontsize=5,
              bbox_to_anchor=(1, 1.1), handletextpad=0.5, columnspacing=0.7)


plt.subplots_adjust(wspace=1.0, hspace=0.6)


#%%