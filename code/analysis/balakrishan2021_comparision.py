#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model 
import diaux.viz 
import diaux.quant
import scipy.stats
cor, pal = diaux.viz.matplotlib_style()

proteome_data = pd.read_csv('../../data/literature/Balakrishnan2021_proteome_fractions.csv')
growth_data = pd.read_csv('../../data/literature/Balakrishnan2021_Fig1A.csv')
growth_data['od600nm'] = np.exp(growth_data['od600nm'])
growth_data['relative_shift_od600nm'] = growth_data['od600nm'] / growth_data[growth_data['relative_shift_time_hr'] >=0]['od600nm'].values[0]
# Infer the growth rates 
post_slc = growth_data.iloc[-4:]
pre_slc = growth_data.iloc[:5]
post_popt = scipy.stats.linregress(post_slc['relative_shift_time_hr'] - post_slc['relative_shift_time_hr'].min(), np.log(post_slc['od600nm']))
pre_popt = scipy.stats.linregress(pre_slc['relative_shift_time_hr'] - pre_slc['relative_shift_time_hr'].min(), np.log(pre_slc['od600nm']))

#%%
# Define the metabolic rates on each substrate using standard parameters to best 
# approximate the growth rates in each regime.
# This was done empirically by fine-tuning nu_max
nu_max = [2.773, 1]

# Run the simulation
suballocation = {'strategy': 'hierarchical',
                 'nu_max': nu_max,
                 'K': [1E-5, 1E-5],
                  'Y': [3E19, 1E19],
                 'n': [1, 1],
                 'hierarchy': [0, 1]}
nutrients = {'init_concs': [0.00061, 0.030]}
species = diaux.model.FluxParityAllocator(suballocation)
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(od_init=0.04, equil_time=100)
species_df, nutrient_df = ecosystem.grow(30, dt=1/60)
steadystates = diaux.quant.profile_steady_states(species_df)
print(steadystates['avg_growth_rate_hr'].values[1])
lagtime, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, steadystates)
lagtime

# Rescale the simulation data
#%%
species_df['relative_shift_mass'] = species_df['M'] / lagtime['stall_biomass'].values[0]
species_df['relative_shift_time_hr'] = species_df['time_hr'] - lagtime['entry_time_hr'].values[0]

#%%
fig, ax = plt.subplots(1,2, figsize=(4, 2))

# ax[0].plot(growth_data['relative_shift_time_hr'], growth_data['relative_shift_od600nm'], 'o',
        # color=cor['primary_black'], ms=4)
ax[0].plot(species_df['relative_shift_time_hr'], np.log(species_df['relative_shift_mass']), '-',
           color=cor['primary_red'], lw=1)

ax[0].set_xlim([-3, 15])
ax[0].set_ylabel('log relative OD$_{600nm}$', fontsize=6)
ax[0].set_ylim(-2, 3)
ax[0].set_xlabel('time from shift [hr]', fontsize=6)
ax[1].plot(species_df['relative_shift_time_hr'], species_df['phi_Mb_1'] + species_df['phi_Mb_2'], '-', color=cor['primary_purple'], lw=1)
ax[1].set_xlim([-3, 15])
_pd = proteome_data[proteome_data['sample'].isin(['5', '15', '60', '120'])]
_pd['sample_time'] = _pd['sample'].astype(float) / 60
ax[1].plot(_pd['sample_time'], _pd['phiMb'], 'o', color=cor['primary_purple'], ms=4)
proteome_data.loc[proteome_data['sample'] == 'preshift', 'sample_time'] = -1
proteome_data.loc[proteome_data['sample'] == 'postshift', 'sample_time'] = 10
proteome_data.dropna(inplace=True)
ax[1].plot(proteome_data['sample_time'], proteome_data['phiMb'], 'o', color=cor['primary_purple'], ms=4)