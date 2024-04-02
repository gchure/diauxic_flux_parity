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
          [0.7, 0.3],
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
feed_concs = [0.0002, 0.0018]
inflow_rates = [dil_rate, dil_rate]
degradation_rates = [dil_rate, dil_rate]
nutrients = {'init_concs': feed_concs,
             'feed_concs': feed_concs,
             'inflow_rates': inflow_rates,
             'degradation_rates': degradation_rates}
ecosystem = diaux.model.Ecosystem(species, nutrients)

#%%
ecosystem.preculture()
species_df, nutrient_df = ecosystem.grow(100, dt=0.001)

#%%
# Assign colors to the species
cmap = sns.husl_palette(n_colors=len(species), l=0.7)
species_df['color'] = [cmap[j-1] for j in species_df['species_label'].values]
# Set up a figure canvas
fig = plt.figure(figsize=(6,2.5))#, layout="")
gs = GridSpec(4, 5)

# Define the gridspec axes
line_ax = fig.add_subplot(gs[0, :2])
bar_ax = fig.add_subplot(gs[1:, :2])
freq_ax = fig.add_subplot(gs[:2, 2:])
nut_ax = fig.add_subplot(gs[2:, 2:])


# Adjust various limits and labels
freq_ax.set_title('ecosystem composition', fontsize=6)
nut_ax.set_title('nutrient composition', fontsize=6)
freq_ax.set_ylabel('mass frequency', fontsize=6)
nut_ax.set_ylabel('concentration [M]', fontsize=6)
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
line_ax.plot([0], [0], 's', markerfacecolor='w', markeredgecolor=cor['primary_black'])
line_ax.text(-0.015, -0.012, 'A', fontsize=4, color=cor['primary_black'])
line_ax.plot([1], [0], 's', markerfacecolor='w', markeredgecolor=cor['primary_black'])
line_ax.text(0.985,  -0.012, 'B', fontsize=4, color=cor['primary_black'])
line_ax.set_title('metabolic strategy 1-simplex', fontsize=6, y=0.8)

bar_ax.set_facecolor('w')
# bar_ax.set_xticks(np.arange(0, len(species), 1))
# bar_ax.set_xticklabels([s.label for s in ecosystem.species])
bar_ax.set_xlabel('species label', fontsize=6)
bar_ax.set_title('metabolic suballocation', fontsize=6, zorder=1000)
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
    line_ax.plot(1-tp['alpha_1'], 0, marker, label=i, ms=4, color=tp['color'])
    bottom = 0
    colors = ['purple', 'light_purple']
    for j in range(2):
        suballoc = d.iloc[-1][f'alpha_{j+1}']
        bar_ax.bar(g, suballoc, bottom=bottom, color=cor[colors[j]])
        bottom += suballoc
    bar_ax.plot(g, -.05, 'o', ms=4, color=tp['color'])
bar_ax.set_xticks(labels)

# Plot the species frequency and nutrient concentrations
for g, d in species_df.groupby(['species_label', 'color']):
    freq_ax.plot(d['time_hr'], d['frequency'], lw=1.5, color=g[1])

# Plot the nutrient breakdown
plt.subplots_adjust(wspace=1.0, hspace=0.5)


#%%