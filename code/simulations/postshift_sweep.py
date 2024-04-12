#%%
import numpy as np 
import pandas as pd
import importlib
import diaux.model
import diaux.quant
import tqdm

# Set the preshift ranges and fixed nutrient concentrations
nu_preshift_range = np.arange(1, 20, 1)
nutrients = {'init_concs': [0.0005, 30]}

# Set the storage dataframes for the trajectories, steadystates, and lagtimes
species_traj = pd.DataFrame([])
steady_states_df =  pd.DataFrame([])
lag_times_df = pd.DataFrame([])
# Set the suballocation
suballocation = {'strategy': 'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'hierarchy': [0, 1]}
# Iterate through each preshift rate
for i, nu_pre in enumerate(tqdm.tqdm(nu_preshift_range, 
                              desc='Iterating through preshift conditions...')):
    # Set and iterate through the postshift rate range
    nu_postshift_range = np.linspace(0.5, nu_pre, 10)
    for j, nu_post in enumerate(tqdm.tqdm(nu_postshift_range,
                             desc='Iterating through postshift conditions...')):
        # Update the preshift/postshift rates
        suballocation['nu_max'] = [nu_pre, nu_post]

        # Define the species and ecosystem
        species = diaux.model.FluxParityAllocator(suballocation)
        ecosystem = diaux.model.Ecosystem([species], nutrients)

        # Preculture to the steady-state regime in the preshift env.
        ecosystem.preculture(verbose=False, init_conc_override=[1.0, 0.0],
                             equil_time=50, max_iter=30, od_init=0.001)

        # Grow the system for the time range with fine time ranges
        species_df, nutrient_df = ecosystem.grow(100, dt=1E-3, verbose=False)

        # compute the steady states
        steady_states = diaux.quant.profile_steady_states(species_df)

        # Profile the lag times
        lag_times, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, 
                                                             species_df,
                                                             steady_states)
        
        # Clip the trajectories to the point where nutrient 2 exhauts.  
        if len(nexh) == 2:
            t_stop = nexh[nexh['nutrient_label']==1]['exhaustion_time_hr'].values[0]
            species_df = species_df[species_df['time_hr'] <= t_stop]

        # Add simulation information to each dataframe
        for d in [species_df, steady_states, lag_times]:
            d['nu_max_preshift'] =  nu_pre
            d['nu_max_postshift'] = nu_post
            d['K'] = suballocation['K'][0]

        # Store the dataframes
        species_traj = pd.concat([species_traj, species_df], sort=False)
        steady_states_df = pd.concat([steady_states_df, steady_states], sort=False)
        lag_times_df = pd.concat([lag_times_df, lag_times], sort=False) 
#%%
species_traj.to_csv('../../data/simulations/preshift_postshift_sweep_species_trajectories.csv', index=False)        
steady_states_df.to_csv('../../data/simulations/preshift_postshift_sweep_steady_states.csv', index=False)
lag_times_df.to_csv('../../data/simulations/preshift_postshift_sweep_lag_times.csv', index=False) 