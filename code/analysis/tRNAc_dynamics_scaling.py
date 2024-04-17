#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import tqdm
import diaux.model 
import diaux.viz 
import diaux.quant
cor, pal = diaux.viz.matplotlib_style()

# Set up storage dataframes for the simulation outputs
df = pd.DataFrame([])
# Define postshift ranges
nu_max_preshift = [1, 3, 5, 7, 10]
n_points = 10

# Define the common suballocation
suballocation = {'strategy': 'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1]}
nutrients = {'init_concs': [0.0005, 10]}

# Iterate through the shift conditions 
for i, nu_pre in enumerate(tqdm.tqdm(nu_max_preshift)):
    nu_max_postshift = np.linspace(0.5, nu_pre, n_points)
    for j, nu_post in enumerate(nu_max_postshift):
        # Instantiate and grow the ecosystem
        suballocation['nu_max'] = [nu_pre, nu_post]
        species = diaux.model.FluxParityAllocator(suballocation)
        ecosystem = diaux.model.Ecosystem([species], nutrients)
        ecosystem.preculture(od_init=0.001, equil_time=100, verbose=False)
        _species_df, _nutrient_df = ecosystem.grow(50, dt=1/60, verbose=False)
        
        # find the time at which nutrient 1 is exhausted
        t_exhaust = _nutrient_df[(_nutrient_df['nutrient_label']==0) & 
                                (_nutrient_df['env_conc'] <= 0)]['time_hr'].values[0]
        M_stall = _species_df[_species_df['time_hr']==t_exhaust]['M'].values[0]

        # Compute the acclimation time.
        postshift = _species_df[_species_df['M'] >= M_stall]
        ind_acc = np.argmax(np.diff(postshift['tRNA_c'].values))
        t_acc = postshift['time_hr'].values[ind_acc+1] - postshift['time_hr'].values[0]

        # Compute the heuristic lag time
        steadystates = diaux.quant.profile_steady_states(_species_df)
        lag_time, _ = diaux.quant.draft_profile_lag_time(_nutrient_df, _species_df, steadystates)
        
        _df = pd.DataFrame({'nu_max_preshift':nu_pre,
                            'nu_max_postshift':nu_post,
                            'lag_time_hr': lag_time['lag_time_hr'].values[0],
                            't_acc':t_acc,
                            'lam_preshift':steadystates['avg_growth_rate_hr'].values[0],
                            'lam_postshift':steadystates['avg_growth_rate_hr'].values[1]},
                            index=[0])
        df = pd.concat([df, _df], sort=False)

