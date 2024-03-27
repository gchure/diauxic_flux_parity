#%%
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import diaux.model
import diaux.quant
import diaux.viz
import tqdm
cor, pal = diaux.viz.matplotlib_style()

nu_post_range = np.linspace(0.5, 10, 20)
trajectories = pd.DataFrame([])
steadystates = pd.DataFrame([])
lagtimes = pd.DataFrame([])
for _, nu_post in enumerate(tqdm.tqdm(nu_post_range)):
    nu_pre_range = np.arange(nu_post + 0.5, 20, 0.5)
    for _, nu_pre in enumerate(nu_pre_range):
        # Define the species which can eat two nutrients
        suballocation = {'strategy': 'dynamic',
                         'nu_max': [nu_pre, nu_post],
                         'hierarchy': [0, 1],
                         'frac_useful': [1, 1],
                         'K': [1E-5, 1E-5],
                         'n': [2, 2]}
        species = diaux.model.FluxParityAllocator(suballocation, metabolic_hierarchy=False,
                                                  label=1)                    
        nutrients = {'init_concs': [0.0005, 100]}
        ecosystem = diaux.model.Ecosystem([species], nutrients)
        ecosystem.preculture(equil_time=100, verbose=False)
        species_df, nutrient_df = ecosystem.grow(100, dt=1/60, verbose=False)
        species_df['nu_max_preshift'] =  nu_pre
        species_df['nu_max_postshift'] = nu_post
        species_df['useful_fraction_pre'] = 1
        species_df['useful_fraction_post'] = 1
        trajectories = pd.concat([trajectories, species_df], sort=False)
        # Profile the steady state
        ss_df = diaux.quant.draft_profile_steady_states(species_df)
        ss_df['nu_max_preshift'] = nu_pre
        ss_df['nu_max_postshift'] = nu_post
        ss_df['useful_fraction_pre'] = 1
        ss_df['useful_fraction_post'] = 1
        steadystates = pd.concat([steadystates, ss_df], sort=False)
        # Compute the lag time
        M_INIT = species_df['M'].values[0]
        M_STAR= (nutrients['init_concs'][0] * species.Y[0] + M_INIT)
        t_hit = species_df.iloc[np.where(species_df['M'] >= M_STAR)[0][0]]['time_hr']
        postshift = ss_df[ss_df['steady_state_idx']==1]
        preshift = ss_df[ss_df['steady_state_idx']==0]
        _lag_phase = species_df[(species_df['time_hr'] >= preshift['time_stop'].values[0]) & 
                          (species_df['time_hr'] <= postshift['time_start'].values[0])]
        lag_time = ((np.log(M_STAR) - np.log(postshift['M_start']))/postshift['growth_rate_hr']) + postshift['time_start']
        _df = pd.DataFrame({'nu_max_preshift':nu_pre,
                            'nu_max_postshift':nu_post,
                            'lag_time_hr': lag_time - t_hit,
                            'shift_time': t_hit,
                            'regrowth_time': lag_time,
                            'shift_biomass': M_STAR,
                            'min_gamma': _lag_phase['gamma'].min(),
                            'preshift_growth_rate_hr':preshift['growth_rate_hr'].values[0],
                            'postshift_growth_rate_hr':postshift['growth_rate_hr'].values[0],
                            'useful_fraction_pre': 1,
                            'useful_fraction_post': 1,
                            },
                            index=[0])
        lagtimes = pd.concat([lagtimes, _df], sort=False) 
#%%
lagtimes.to_csv('../../data/simulations/preshift_variation_lagtimes.csv', index=False)
steadystates.to_csv('../../data/simulations/preshift_variation_steadystates.csv', index=False)
trajectories.to_csv('../../data/simulations/preshift_variation_trajectories.csv', index=False)
