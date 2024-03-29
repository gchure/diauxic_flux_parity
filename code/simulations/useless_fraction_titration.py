#%%
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import diaux.model
import diaux.quant
import tqdm

nu_pre = 10.0
nu_post_range = np.linspace(1, 9, 10)
step = 0.01
useful_fracs = np.linspace(0.1, 1, 10)
trajectories = pd.DataFrame([])
steadystates = pd.DataFrame([])
lagtimes = pd.DataFrame([])

#%%
# Define the species which can eat two nutrients
for i, f in enumerate(tqdm.tqdm(useful_fracs)):
    for j, nu_post in enumerate(nu_post_range):
        suballocation = {'strategy': 'dynamic',
                         'nu_max': [nu_pre, nu_post],
                         'frac_useful': [1, f],
                         'K': [1E-5, 1E-5],
                         'n': [1, 1]}
        species = diaux.model.FluxParityAllocator(suballocation, metabolic_hierarchy=False, label=1)                    
        nutrients = {'init_concs': [0.001, 100]}
        ecosystem = diaux.model.Ecosystem([species], nutrients)
        ecosystem.preculture(equil_time=100, verbose=False)
        species_df, nutrient_df = ecosystem.grow(500, dt=1/60, verbose=False)
        species_df['nu_max_preshift'] =  nu_pre
        species_df['nu_max_postshift'] = nu_post
        species_df['useful_fraction_pre'] = 1
        species_df['useful_fraction_post'] = f
        trajectories = pd.concat([trajectories, species_df], sort=False)

        # Profile the steady state
        ss_df = diaux.quant.draft_profile_steady_states(species_df)
        ss_df['nu_max_preshift'] = nu_pre
        ss_df['nu_max_postshift'] = nu_post
        ss_df['useful_fraction_pre'] = 1
        ss_df['useful_fraction_post'] = f
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
                                'useful_fraction_post': f,
                                },
                                index=[0])
        lagtimes = pd.concat([lagtimes, _df], sort=False) 

#%%
lagtimes.to_csv('../../data/simulations/useful_frac_sweep_lagtimes.csv', index=False)
steadystates.to_csv('../../data/simulations/useful_frac_sweep_steadystates.csv', index=False)
trajectories.to_csv('../../data/simulations/useful_frac_sweep_trajectories.csv', index=False)