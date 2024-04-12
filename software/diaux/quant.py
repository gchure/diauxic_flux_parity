#%%
import numpy as np 
import pandas as pd
from . import model


def profile_steady_states(species_df, 
                     alloc_stability_thresh=5E-4,
                     tRNA_stability_thresh=1E-3,
                     growth_stability_thresh=5E-5):
    """
    Identifies the steady-state phases of growth observed in each dilution cycle 
    of a simulation for each individual species.

    Parameters
    ----------
    species_df : pandas DataFrame 
        A Pandas DataFrame containing the cell physiological details of the 
        simulated growth.
    alloc_stability_thresh: float
        The threshold below which the ratio of the ribosomal allocation and 
        the ribosome content are considered to be the same. Default is 0.005.
    tRNA_stability_thresh: float 
        The threshold by which the difference in the charged-tRNA concentration 
        stability  is considered to be stationary. Default is 0.005.

    Returns
    -------
    steady_states : pandas DataFrame   
        A pandas DataFrame containing information for the start, stop, duration, 
        and other statistics about the steadystates shifts observed in the simulation 
        for each species and each dilution cycle.
    """
    steady_states = pd.DataFrame([])

    # Select the steadystates
    # _species_df = species_df.copy(deep=True)
    # frac_alloc_ss = np.abs(species_df['alloc_stability'])
    # alloc_ss = np.abs(np.diff(1 - species_df['alloc_stability'])) <= alloc_stability_thresh
    # tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= tRNA_stability_thresh
    # _species_df['steady_state'] = alloc_ss * tRNA_c_ss

    # Iterate through each species and dilution phase
    for g, d in species_df.groupby('species_label'):
        n_ss_total = 0
        for _g, _d in d.groupby('dilution_cycle'):
            lam = _d['gamma'] * _d['ribosome_content']
            frac_dlam = np.abs(np.diff(lam)/lam[:-1])
            frac_dalloc = np.abs(np.diff(_d['alloc_stability'].values)/_d['alloc_stability'].values[:-1])
            _d['steady_state'] = (frac_dlam <= growth_stability_thresh)# * (frac_dalloc <= alloc_stability_thresh) 
            # Find the indices where switchpoints happen, indicating a transition
            # between states
            ss_vals = _d['steady_state'].values
            switch_inds = np.where(np.diff(ss_vals) != False)[0]
            if len(switch_inds) > 0:
                # Set up storage lists for finding the start and stop points for 
                # each phase
                start_stop = [[], []]
                for i, idx in enumerate(switch_inds):
                    # Determine if the switch is in or out of steady-state
                    if ss_vals[idx+1] == False:
                        # If this is the first ind, this means that the simulation
                        # at this phase *starts* in steady state, so keep track 
                        # of the idx as 0
                        if i == 0:
                            start_stop[0].append(0) 
                        start_stop[1].append(idx)  
                    # If not False, then by definition a steady-state is beginning
                    else:
                        start_stop[0].append(idx)
                # Determine if there are more starts then stops. If that's the case,
                # the steady state must end at the maximum idx, meaning the state 
                # only ends because the simulation ends here. 
                if len(start_stop[0]) > len(start_stop[1]):
                    start_stop[1].append(-1)

                # Iterate through each phase, given the start and stop indexes,
                # isolate the, phase, and compute the interesting properties 
                for i, (start_idx, stop_idx) in enumerate(zip(*start_stop)):
                    n_ss_total += 1
                    _start = _d.iloc[start_idx]
                    _stop = _d.iloc[stop_idx]

                    ss_phase = _d.iloc[start_idx:stop_idx]
                    ss_lam = ss_phase['gamma'] * ss_phase['ribosome_content']
                    _df = pd.DataFrame({'time_start': _start['time_hr'],
                                        'time_stop': _stop['time_hr'],
                                        'M_start':_start['M'],
                                        'M_stop':_stop['M'],
                                        'avg_growth_rate_hr': np.mean(ss_lam),
                                        'min_growth_rate_hr': np.min(ss_lam),
                                        'max_growth_rate_hr' : np.max(ss_lam),
                                        'dilution_cycle':_g,
                                        'species_label':g,
                                        'cycle_ss_idx': i+1,
                                        'total_ss_idx': n_ss_total},
                                    index=[0])
                    steady_states = pd.concat([steady_states, _df], sort=False)
    return steady_states


def draft_profile_lag_time(nutrient_df,
                           species_df,
                           steady_states):
    """
    Computes the lag time between a preshift and postshift nutrient pair. 
    TODO: Write this in a more general way which used information from the 
    hierarchy to figure out the pairing 
    """

    # Find the shift biomasses
    nutrient_exhaustion = pd.DataFrame([])
    for g, d in nutrient_df.groupby('dilution_cycle'):
        for _g, _d in d.groupby('nutrient_label'):
            exhaust_ind = np.where(_d['env_conc']<=0)[0]
            if len(exhaust_ind) > 0:
                exhaust_time = _d.iloc[exhaust_ind[0]]['time_hr']
                for __g, __d in species_df[species_df['dilution_cycle']==g].groupby('species_label'):
                    exhaust_biomass = __d.iloc[exhaust_ind[0]]['M']
                    _df = pd.DataFrame({'dilution_cycle': g,
                                        'nutrient_label': _g,
                                        'species_label': __g,
                                        'exhaustion_time_hr': exhaust_time,
                                        'exhaustion_biomass': exhaust_biomass},
                                        index=[0])
                    nutrient_exhaustion = pd.concat([nutrient_exhaustion, _df], sort=False)
    # Find the lag times
    lag_times = pd.DataFrame([])
    for g, d in nutrient_exhaustion.groupby(['dilution_cycle', 'species_label', 'nutrient_label']):
        # Find the species steadystate pair where the exhaustion time is 
        # interleaved 
        t_exh = d['exhaustion_time_hr'].values[0]
        M_stall = d['exhaustion_biomass'].values[0]
        ss_phase = steady_states[(steady_states['dilution_cycle']==g[0]) & 
                                 (steady_states['species_label']==g[1]) &
                                 (steady_states['time_start'] >= t_exh)]
        if len(ss_phase) > 0:
           # Compute the heuristic for the lag time
           log_mass_diff = np.log(M_stall) - np.log(ss_phase['M_start'].values[0])
           postshift_lam = ss_phase['avg_growth_rate_hr'].values[0]
           t_star = ss_phase['time_start'] + log_mass_diff/postshift_lam 
           lag_time = t_star - t_exh
    
           # Assemble the dataframe
           _df = pd.DataFrame({'dilution_cycle':g[0],
                               'species_label':g[1],
                               'preshift_nutrient_label':g[2],
                               'stall_biomass':M_stall,
                               'entry_time_hr': t_exh,
                               'exit_time_hr': t_star,
                               'lag_time_hr':lag_time})
           lag_times = pd.concat([lag_times, _df], sort=False)
    return [lag_times, nutrient_exhaustion]



def draft_profile_steady_states(species_df, 
                                alloc_stability_thresh=5E-3,
                                tRNA_stability_thresh=1E-3):
    """
    DEPRECATED. Use `profile_steady_states` instead 
    """
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
