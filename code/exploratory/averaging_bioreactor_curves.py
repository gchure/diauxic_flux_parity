#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import diaux.model 
import diaux.viz
import diaux.quant
cor, pal = diaux.viz.matplotlib_style()

data =  pd.read_csv('../../data/bioreactor/pilot_diauxic_shifts.csv')
lag_data = pd.read_csv('../../data/bioreactor/pilot_lag_times.csv')
avg_preshift = lag_data['preshift_growth_rate'].mean()
avg_postshift = lag_data['postshift_growth_rate'].mean()
print(avg_preshift, avg_postshift)
data['inst_growth_rate'] = data['log_od650nm_filtered'].diff() / data['time_hr'].diff()

#%% 

# Set up the model
suballocation = {'strategy': 'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'Y': [3E19, 1E19],
                 'nu_max': [3.2, 1.5]}
                 
nutrients = {'init_concs': [0.0005, 0.030]}
species = diaux.model.FluxParityAllocator(suballocation)
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(equil_time=100)
species_df, nutrient_df = ecosystem.grow(30, dt=1/60)
steadystates = diaux.quant.profile_steady_states(species_df)
lag_times, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, steadystates)
species_df['relative_shift_mass'] = species_df['M'] / lag_times['stall_biomass'].values[0]
species_df['relative_shift_time_hr'] = species_df['time_hr'] - lag_times['entry_time_hr'].values[0]
species_df['inst_growth_rate'] = species_df['gamma'] * species_df['ribosome_content']
print(steadystates['avg_growth_rate_hr'].values)
# Convert time from hours to minutes
data['time_min_centered'] = data['time_hr_centered'] * 60

# Create time windows of 5 minutes
data['time_window'] = np.floor(data['time_min_centered'] / 5) * 5

data = data[data['time_hr_centered'] > -1]

# Group by time window and calculate the mean of the signal
averaged_data = data.groupby('time_window')[['od650nm_centered', 'inst_growth_rate']].median().reset_index()
averaged_data['time_hr_centered'] = averaged_data['time_window'] / 60
fig, ax = plt.subplots(2,1, figsize=(4, 4), sharex=True)
ax[0].plot(data['time_hr_centered'], np.log(data['od650nm_centered']), 'k.',
         markeredgewidth=0, ms=3, alpha=0.2)
ax[0].plot(averaged_data['time_hr_centered'], np.log(averaged_data['od650nm_centered']), 
          color=cor['pale_blue'], marker='o', ms=3, linestyle='none',
          markeredgewidth=0.5, markeredgecolor=cor['blue'])
ax[0].plot(species_df['relative_shift_time_hr'], np.log(species_df['relative_shift_mass']), '-',
            color=cor['primary_red'], lw=1)
ax[1].plot(data['time_hr_centered'], (data['inst_growth_rate']), 'k.',
         markeredgewidth=0, ms=3, alpha=0.2)
ax[1].plot(averaged_data['time_hr_centered'], averaged_data['inst_growth_rate'], 
          color=cor['pale_blue'], marker='o',markeredgewidth=0.5, 
          markeredgecolor=cor['blue'], ms=3, linestyle='none')
ax[1].plot(species_df['relative_shift_time_hr'], species_df['inst_growth_rate'], '-',
           color=cor['primary_red'], lw=1)
ax[1].set_ylim([-0.3, 1.2])
ax[0].set_xlim([-1.2, 4.5])
ax[0].set_ylim([-1, 1])