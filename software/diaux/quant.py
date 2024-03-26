#%%
import numpy as np 
import pandas as pd
from . import model

def draft_profile_steady_states(species_df, 
                                alloc_stability_thresh=5E-3,
                                tRNA_stability_thresh=1E-3):
    # Identify the steady-state regimes
    alloc_stability_thresh = 5E-3
    tRNA_stability_thresh = 1E-3
    alloc_ss = np.abs(1 - species_df['alloc_stability']) <= alloc_stability_thresh
    tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= tRNA_stability_thresh
    species_df['steady_state'] = alloc_ss * tRNA_c_ss
    ss_df = pd.DataFrame([])
    for g, d in species_df.groupby('species_label'):
        ss_v = d['steady_state'].values
        _diff = np.diff(ss_v)
        inds = np.where(_diff != False)[0]
        if len(inds) > 0:
            start, stop = [], []
            for i, idx in enumerate(inds):  
                if ss_v[idx+1] == False:
                    if (i==0):
                      start.append(0)
                    stop.append(idx)
                else:
                    start.append(idx)
            if len(start) > len(stop):
                stop.append(len(ss_v))

            for i, (_start, _stop) in enumerate(zip(start, stop)):
                __start = species_df.iloc[_start]
                __stop = species_df.iloc[_stop]
                phase = d.iloc[_start:_stop]
                phase_lam = phase['gamma'].values * phase['ribosome_content'].values
                _df = pd.DataFrame({'time_start':__start['time_hr'],
                                    'time_stop':__stop['time_hr'],
                                    'M_start':__start['M'],
                                    'M_stop':__stop['M'],
                                    'growth_rate_hr': phase_lam.mean(),
                                    'species_label': g,
                                    'steady_state_idx': i},
                                    index=[0]) 
                ss_df = pd.concat([_df, ss_df], sort=False)
    return ss_df 
