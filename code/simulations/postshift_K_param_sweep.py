#%%
import numpy as np 
import pandas as pd
import diaux.model
import diaux.quant
import diaux.viz
import tqdm
cor, pal = diaux.viz.matplotlib_style()
K_range = np.logspace(-9, 0, 10)
nu_preshift = 10
nu_postshift = np.linspace(1, 20, 20)
nutrients = {'init_concs':[0.0005, 10]}

lag_df = pd.DataFrame([])
for i, K in enumerate(tqdm.tqdm(K_range)):
    for j, nu_post in enumerate(nu_postshift):
        suballocation = {'strategy': 'hierarchical',
                         'K': [K, 1E-5],
                         'n': [1, 1],
                         'nu_max':[nu_preshift, nu_post]}
        species = diaux.model.FluxParityAllocator(suballocation)
        eco = diaux.model.Ecosystem([species], nutrients)
        eco.preculture(verbose=False, init_conc_override=[1, 0])
        species_df, nutrient_df = eco.grow(50, dt=1E-3, verbose=False)
        steadystates = diaux.quant.profile_steady_states(species_df)
        steadystates['duration'] = steadystates['time_stop'] - steadystates['time_start']
        steadystates = steadystates[steadystates['duration']>=0.1]
        lagtimes, _ = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, steadystates)
        lagtimes['nu_preshift'] = nu_preshift
        lagtimes['nu_postshift'] = nu_post
        lagtimes['K'] = K
        lagtimes['lam_postshift'] = steadystates['avg_growth_rate_hr'].values[-1]
        lagtimes['lam_preshift'] = steadystates['avg_growth_rate_hr'].values[0]
        lagtimes['alpha_2_preallocation'] = species_df[(species_df['time_hr']>= steadystates.iloc[0]['time_start']) & 
                                                       (species_df['time_hr']<=steadystates.iloc[0]['time_stop'])]['alpha_2'].mean()
        lag_df = pd.concat([lag_df, lagtimes], sort=False)
#%%
lag_df.to_csv('../../data/simulations/postshift_K_param_sweep_lagtimes.csv')