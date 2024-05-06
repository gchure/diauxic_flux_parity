"""
This script is used to compare the lag time and acclimation time of lag time 
simulations across a breadth of pre and post shifts.
"""
#%%
import tqdm
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.model
import diaux.quant
import diaux.viz 
cor, pal = diaux.viz.matplotlib_style()

# Define preshift ranges
nu_max_preshift = np.arange(1, 10, 1)

# Set up the suballocation
suballocation = {'strategy': 'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'Y': [3E19, 3E19],
                 'hierarchy': [0, 1]}
nutrients = {'init_concs': [0.0005, 1.0]}
acc_df = pd.DataFrame([])
for j, nu_pre in enumerate(tqdm.tqdm(nu_max_preshift)):
    nu_max_postshift = np.arange(0.5, nu_pre+0.5, 0.05)
    for j, nu_post in enumerate(tqdm.tqdm(nu_max_postshift)):
        suballocation['nu_max'] = [nu_pre, nu_post]
        species = diaux.model.FluxParityAllocator(suballocation)
        ecosystem = diaux.model.Ecosystem([species], nutrients)
        ecosystem.preculture(equil_time=100, verbose=False)
        species_df, nutrient_df = ecosystem.grow(50, dt=1/60, verbose=False)
        steadystates = diaux.quant.profile_steady_states(species_df)
        lag_times, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, steadystates)

        # Compute the tRNA acclimation time
        tRNA_t_acc = diaux.quant.compute_acclimation_time(species_df, nutrient_df)
        species_df['lam'] = species_df['gamma'] * species_df['ribosome_content']
        lam_t_acc = diaux.quant.compute_acclimation_time(species_df, nutrient_df,
                                                          quantity='lam')
        gamma_t_acc = diaux.quant.compute_acclimation_time(species_df, nutrient_df,
                                                          quantity='gamma')

        _acc_df = pd.DataFrame({'nu_max_preshift': nu_pre,
                               'nu_max_postshift': nu_post,
                               'lag_time': lag_times['lag_time_hr'].values[0],
                               'tRNAc_acclimation_time': tRNA_t_acc,
                               'lambda_acclimation_time': lam_t_acc,
                               'gamma_acclimation_time': gamma_t_acc,
                               'preshift_growth_rate_hr': steadystates['avg_growth_rate_hr'].values[0],
                               'postshift_growth_rate_hr': steadystates['avg_growth_rate_hr'].values[1],},
                               index=[0])
        acc_df = pd.concat([acc_df, _acc_df], sort=False)
#%%
acc_df.to_csv('../../data/simulations/acclimation_times.csv', index=False)
