#%%
import numpy as np 
import pandas as pd 
import diaux.model
import diaux.quant
import tqdm

nu_pre = 5.0
nu_post = 1.0
alpha_range = np.linspace(0.1, 1, 10)
nutrients = {'init_concs': [0.001, 100]}
trajectories = pd.DataFrame([])
steadystates = pd.DataFrame([])
for i, alpha in enumerate(tqdm.tqdm(alpha_range)):
    suballocation = {'strategy': 'static',
                     'alpha':[alpha, 1-alpha],
                     'nu_max': [nu_pre, nu_post],
                     'frac_useful': [1, 1],
                     'K': [1E-5, 1E-5],
                     'n': [1, 1]}
    species = diaux.model.FluxParityAllocator(suballocation, label=1)                    
    ecosystem = diaux.model.Ecosystem([species], nutrients)
    ecosystem.preculture(equil_time=300, verbose=False)
    species_df, nutrient_df = ecosystem.grow(100, dt=1/60, verbose=False) 
    species_df['alpha_1'] = alpha 
    species_df['alpha_2'] = 1 - alpha
    species_df['nu_max_preshift'] = nu_pre 
    species_df['nu_max_postshift'] = nu_post
    trajectories = pd.concat([trajectories, species_df], sort=False)

    ss_df = diaux.quant.draft_profile_steady_states(species_df)
    ss_df['alpha_1'] = alpha
    ss_df['alpha_2'] = 1 - alpha
    ss_df['nu_max_preshift'] = nu_pre
    ss_df['nu_max_postshift'] = nu_post
    steadystates = pd.concat([steadystates, ss_df], sort=False)

#%%
trajectories.to_csv('../../data/simulations/static_suballocation_trajectories.csv', index=False)
steadystates.to_csv('../../data/simulations/static_suballocation_steadystates.csv', index=False)