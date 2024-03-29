#%%
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import diaux.model
import diaux.quant
import tqdm

nu_pre = 10.0
nu_post = 1.0


suballocation = {'strategy': 'dynamic',
                     'nu_max': [nu_pre, nu_post],
                     'frac_useful': [1, 1],
                     'K': [1E-5, 1E-5],
                     'n': [1, 1]}
species = diaux.model.FluxParityAllocator(suballocation, metabolic_hierarchy=False, label=1)                    
nutrients = {'init_concs': [0.001, 100]}
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(equil_time=100, verbose=False)
species_df, nutrient_df = ecosystem.grow(500, 
                                         dt=1/60, 
                                         bottleneck={'type':'time',
                                                     'interval':30,
                                                     'target':0.04})
#%%
steadystates = pd.DataFrame([])
lagtimes = pd.DataFrame([])
for g, d in tqdm.tqdm(species_df.groupby('dilution_cycle')): 
    # Profile the steady state
    d['time_hr'] -= d['time_hr'].min()
    ss_df = diaux.quant.draft_profile_steady_states(d)
    ss_df['nu_max_preshift'] = nu_pre
    ss_df['nu_max_postshift'] = nu_post
    ss_df['useful_fraction_pre'] = 1
    ss_df['useful_fraction_post'] = 1
    ss_df['dilution_cycle'] = g
    steadystates = pd.concat([steadystates, ss_df], sort=False)

    # Compute the lag time
    M_INIT = d['M'].values[0]
    M_STAR= (nutrients['init_concs'][0] * species.Y[0] + M_INIT)
    t_hit = d.iloc[np.where(d['M'] >= M_STAR)[0][0]]['time_hr']

    if len(ss_df) > 1:
        ss_df = ss_df[ss_df['steady_state_idx']==1]
    _lag_phase = d[(d['M'] >= M_STAR)]
                      
    lag_time = ((np.log(M_STAR) - np.log(ss_df['M_start']))/ss_df['growth_rate_hr']) + ss_df['time_start']
    _df = pd.DataFrame({'nu_max_preshift':nu_pre,
                        'nu_max_postshift':nu_post,
                        'lag_time_hr': lag_time,
                        'min_gamma': _lag_phase['gamma'].min(),
                        'useful_fraction_pre': 1,
                        'useful_fraction_post': 1,
                        'dilution_cyle':g},
                        index=[0])
    lagtimes = pd.concat([lagtimes, _df], sort=False) 

#%%
lagtimes.to_csv('../../data/simulations/preallocation_lagtimes.csv', index=False)
