#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model 
import diaux.viz 
import diaux.quant
import importlib
importlib.reload(diaux.quant)

cor, pal = diaux.viz.matplotlib_style()

# Define a representative growth curve with annotated lag times

suballocation = {'strategy': 'hierarchical',
                 'nu_max': [10.0, 1.0],
                 'n': [1, 1],
                 'K':[1E-5, 1E-5]}
nutrients = [0.0005, 0.030]
od_init = 0.001
species = diaux.model.FluxParityAllocator(suballocation)
ecosystem = diaux.model.Ecosystem([species], {'init_concs':nutrients})
ecosystem.preculture(od_init=od_init)
stall_biomass = od_init * diaux.model.OD_CONV + nutrients[0] * species.Y[0]
species_df, nutrient_df = ecosystem.grow(20, dt=1E-3)

#%%
steadystates = diaux.quant.profile_steady_states(species_df)
steadystates

#%%
# Compute the lag times
lag_times, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, steadystates)

#%%

# Instantiate the figure canvas and set appropriate labeling
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_ylabel('OD$_{600nm}$ / mL\napproximate optical density', fontsize=6)
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_yscale('log')
ax.set_xlim([0, 20])
ylims = [0.001, 5.5]

# Add shaded regions highlighting the two steady-state regimes
ax.plot(species_df['time_hr'], species_df['M']/diaux.model.OD_CONV, 'k-',
         lw=1, zorder=100)

# Shade in the steadystate time ranges
for g, d in steadystates.groupby('total_ss_idx'):
    ax.fill_betweenx(ylims, d['time_start'], d['time_stop'], color=cor['pale_black'],
                     alpha=0.35)
    # Add labels
    if g == 1:
        y = 3
        tshift = 0
        txt = 'SS$_{pre}$'
    else:
        y = 0.0015
        txt = 'SS$_{post}$'
    ax.text(d['time_start'].values[0]+0.35, y, txt, fontsize=5, 
                color=cor['light_black'])

# Add a line indicating the biomass point where the growth stalled
ax.hlines(stall_biomass/diaux.model.OD_CONV, species_df['time_hr'].values[0],
          species_df['time_hr'].values[-1], color=cor['primary_blue'], linestyle='--',
          lw=1, zorder=101)

ax.text(13.5, 0.08, 'stall biomass', fontsize=5, 
        color=cor['primary_blue'])

# Make a plot of the post-shift exponential growth rate
time_range = np.linspace(0, 15, 100)

fit = stall_biomass*np.exp(steadystates['avg_growth_rate_hr'].values[1] * time_range)
ax.plot(time_range + lag_times['exit_time_hr'].values[0],
        fit/diaux.model.OD_CONV, '--', color=cor['primary_green'], lw=1, zorder=110)
ax.text(11, 0.35, 'postshift exponential growth', fontsize=5, color=cor['primary_green'],
        rotation=28,zorder=1000)

# Shade and label the lag time
ax.text(lag_times['entry_time_hr'].values[0]+0.35, 3, 'lag phase', fontsize=5,
        color=cor['primary_blue'])
ax.fill_betweenx(ylims, lag_times['entry_time_hr'], lag_times['exit_time_hr'],
                 color=cor['light_blue'], alpha=0.5)

ax.set_ylim(ylims)
plt.savefig('../../figures/doc/annotated_lag_time.pdf')