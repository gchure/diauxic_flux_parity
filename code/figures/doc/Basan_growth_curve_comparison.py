#%%
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 
import diaux.quant
import diaux.viz 
import diaux.model
cor, pal = diaux.viz.matplotlib_style()
data = pd.read_csv('../../../data/literature/Basan2020_figx1b.csv')
# Set the simulation parameters
suballocation = {'strategy':'hierarchical',
                 'nu_max': [2, 1.5],
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'Y': [3E19, 1E19],
                 'hierarchy': [0, 1]}
nutrients = {'init_concs':[0.0017, 0.060]}
species = diaux.model.FluxParityAllocator(suballocation)
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(od_init=0.001, equil_time=100)
species_df, nutrient_df = ecosystem.grow(100, dt=1/60)
steadystates = diaux.quant.profile_steady_states(species_df)
lagtimes, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, 
                                                    steadystates)

fig, ax = plt.subplots(1,1, figsize=(3, 2))
ax.set_title('Basan et al. 2020: glucose-acetate diauxic shift', fontsize=6)
for g, d in data.groupby('method'):
    if g == 'diauxie':
        m = 'o'
    else:
        m = 's' 
    ax.plot(d['relative_time_hr'], np.log(d['relative_od']), marker=m, ms=5, linestyle='none',
            alpha=0.95, label=f'{g} growth protocol')

# Rescale the simulation
species_df['OD'] = species_df['M']/diaux.model.OD_CONV
species_df['relative_biomass'] = species_df['M'] / lagtimes['stall_biomass'].values[0]
species_df['relative_time'] = species_df['time_hr'] - lagtimes['entry_time_hr'].values[0]

ax.plot(species_df['relative_time'], np.log(species_df['relative_biomass']), '-', 
        lw=1, color=cor['primary_red'], label='prediction')
ax.legend()
ax.set_ylim(-1.5, 3)
ax.set_xlim(-3, 10)
ax.set_ylabel('log relative optical density', fontsize=6)
ax.set_xlabel('time from shift [hr]', fontsize=6)
species_df.to_csv('./basan_simulation.csv')
# #%%
# import scipy.stats
# preshift_data = data[data['relative_time_hr'] <= 0]
# postshift_data = data[data['relative_time_hr'] >= 6]
# preshift_lam = scipy.stats.linregress(preshift_data['relative_time_hr']-preshift_data['relative_time_hr'].values[0],
#                                        np.log(preshift_data['relative_od']))[0]

# postshift_lam = scipy.stats.linregress(postshift_data['relative_time_hr']-postshift_data['relative_time_hr'].values[0],
#                                        np.log(postshift_data['relative_od']))[0]

# postshift_lam