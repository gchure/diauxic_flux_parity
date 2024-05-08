#%%
import numpy as np 
import pandas as pd 
from . import model
from scipy.signal import savgol_filter as scipy_savgol
from scipy.stats import pearsonr as scipy_pearsonr
from scipy.stats import linregress as scipy_linregress

def compute_pearson_correlation(data, 
                                time_col='time_hr',
                                od_col='od650nm',
                                pearson_window=30,
                                savgol_filter=True,
                                savgol_window_length=30,
                                savgol_polyorder=1):
    """
    Given a growth curve compute a rolling Pearson correlation on a filtered 
    signal. 

    Parameters
    ----------
    data : pandas.DataFrame
        The growth curve data that should be processed.
    time_col : str
        The column name of the time dimension in the DataFrame. Default is `time_hr.`
    od_col : str
        The column name of the optical density dimension in the DataFrame. Default 
        is `od650nm`.
    pearson_window: integer
        The size of the window over which to compute the correlation coefficient. 
        Default is 30 points
    savgol_filter : bool
        If True, the growth curve will be smoothed with a Savitsky-Golay filter
        before the correlation coefficient is applied. Default is True 
    savgol_window_length: integer
        The size of the kernel for the Savitsky-Golay filter. Only applied if
        `savgol_filter` is True. Default is 50 time points. 
     savgol_polyorder: integer
        The polynomial order of the Savitsky-Golay filter. Default is 1. Only 
        applied if `savgol_filter` is True.

    Returns
    -------
    corr_df : pandas DataFrame
        A pandas DataFrame with the correlation for each window along with the 
        filtered log optical density, if the Savitsky-Golay filter is applied.
    """
    # Unpack the time vector
    time = data[time_col].values

    # Log transform the the optical density and median filt if desired. 
    log_od = np.log(data[od_col].values)

    if savgol_filter:
        log_od = scipy_savgol(log_od, window_length=savgol_window_length,
                               polyorder=savgol_polyorder)
        log_label = '_filtered'
    else:
        log_label = ''

    # Compute the pearson correlation coefficients and apply the threshold     
    corr = np.empty(len(log_od) - pearson_window)
    for i in range(len(log_od)-pearson_window):
        corr[i] = scipy_pearsonr(time[i:i+pearson_window], log_od[i:i+pearson_window])[0]

    # Copy and restrict the dataframe
    corr_df = data[:-pearson_window].copy()
    corr_df[f'log_{od_col}{log_label}'] = log_od[:-pearson_window]
    corr_df['pearson_correlation'] = corr                                       
    return corr_df


def classify_diauxic_phases(data, 
                          time_col='time_hr',
                          corr_col='pearson_correlation',
                          od_col = 'od650nm',
                          center = True,
                          exponential_pearson_thresh=0.99,
                          stall_pearson_thresh=0,
                          trim=True,
                          minimum_exponential_length=30,
                          ):
    """
    Given a growth curve with computed pearson correlations, annotate the 
    different phases.

    Parameters
    ----------
    data : pandas.DataFrame
        The data with the computed pearson correlation over the growth phases.
    time_col : str
        The column name of the time dimension in the DataFrame. Default is `time_hr.`
    corr_col : str
        The column name of the pearson correlation coefficient in the DataFrame.  
        Default is `pearson_correlation.`
    od_col : str
        The column name of the optical density in the DataFrame.  This is only 
        used if `center` is `True. 
    center : bool
        If True, the DataFrame will be centered in time and optical density to
        the onset of the shift.
    exponential_pearson_thresh: float, [-1, 1]
        The threshold above which the Pearson correlation of the log optical 
        density is considered to be linear, and therefore in expoenential
        growth. 
    stall_pearson_thresh: float
        The threshold below which the Pearson correlation coefficient of the 
        log optical density is considered to be in a stall phase.
    trim : bool
        If True, the data will be trimmed to begin at the beginning of the 
        preshift exponential phase and end at the final point of the postshift
        exponential phase.
    minimum_exponential_length : int
        The minimum number of contiguous timepoints that a Pearson correlation 
        above the `exponential_pearson_thresh` can be considered to be a true 
        exponential phase growth. Default is 10 time points. 

    Returns
    -------
    labeled_df: pandas.DataFrame
       The supplied data with the various phases labeled. Labels include `preshift_exponential`,
       `postshift_exponential`, `stall`, and `transition`.  
    """
    data = data.copy(deep=True)
    # Unpack the correlation to a vector to avoid pandas slowness
    corr = data[corr_col].values
    time = data[time_col].values

    # Apply the thresholds
    stall = corr <= stall_pearson_thresh
    exp = corr >= exponential_pearson_thresh   
    _exp_switchpoints = np.where(np.diff(np.sign(exp.astype(float))) != 0)[0]
    stall_switchpoints = np.where(np.diff(np.sign(stall.astype(float))) != 0)[0]

    # Pad the exponential switchpoints
    if _exp_switchpoints[0] != 0:
        _exp_switchpoints = np.insert(_exp_switchpoints, 0, 0)
    _exp_switchpoints = np.append(_exp_switchpoints, len(corr)-1)     

    # Compute the phase lengths and fold switchpoints that don't meet the 
    # criteria
    lengths = np.diff(_exp_switchpoints)
    exp_switchpoints = [_exp_switchpoints[0]]
    for i, ell in enumerate(lengths):
        if ell > minimum_exponential_length:
            exp_switchpoints.append(_exp_switchpoints[i+1])

    # Assign the default phases. 
    data['phase'] = 'transition'

    # Assign the stall phase    
    stall_start = time[stall_switchpoints[0]-1]
    stall_stop = time[stall_switchpoints[-1]-1]

    data.loc[(data[time_col] >= stall_start) & 
             (data[time_col] <= stall_stop), 'phase'] = 'stall'

    intervals = {}
    if exp[0]:
        intervals['preshift']  = [time[exp_switchpoints[0]], time[exp_switchpoints[1]]] 
    else:
        intervals['preshift']  = [time[exp_switchpoints[1]], time[exp_switchpoints[2]]] 
    if exp[1]:
        intervals['postshift'] = [time[exp_switchpoints[-2]] + 1, time[exp_switchpoints][-1]]
    else: 
        intervals['postshift'] = [time[exp_switchpoints[-2]] + 1, time[exp_switchpoints][-1]]
    for k, v in intervals.items():
        data.loc[(time >= v[0]) & (time <= v[1]), 'phase'] = f'{k}_exponential'

    if trim:
        data = data[(data[time_col] >= intervals['preshift'][0]) & 
                    (data[time_col] <= intervals['postshift'][1])]

    # Add centered columns if desired
    if center:
        stall_time_init = data[data['phase']=='stall'][time_col].values[0]
        stall_od_init = data[data['phase']=='stall'][od_col].values[0]
        data[time_col +'_centered'] = data[time_col] - stall_time_init
        data[od_col + '_centered'] = data[od_col].values / stall_od_init 
        if f'log_{od_col}_filtered' in data.keys():
            data[f'log_{od_col}_filtered_centered'] = data[f'log_{od_col}_filtered'].values / data[data[time_col + '_centered']==0][f'log_{od_col}_filtered'].values[0]
    return data

def quantify_lag_time(data,
                      time_col='time_hr',
                      od_col='od650nm'):
    """
    Calculates the heuristic lag time from a labeled diauxic growth curve.

    Parameters
    ----------
    data : pandas.DataFrame
        A Pandas DataFrame containing the processed and labeled growth curve.
    time_col: str
        The column name of the time dimension in the DataFrame. Default is `time_hr.`
    od_col: str
        The column name of the optical density in the DataFrame. Default is `od650nm.`

    Returns
    -------
    lag_time : pandas.DataFrame
        A pandas DataFrame reporting the lag time, along with the average stall 
        biomass and the preshift/postshift exponential growth rates.
    """

    # Fit the exponential growth rates for each exponential phase
    preshift = data[data['phase']=='preshift_exponential']
    postshift = data[data['phase']=='postshift_exponential']
    preshift_popt = scipy_linregress(preshift[time_col]-preshift[time_col].min(),
                                   np.log(preshift[od_col]))
    postshift_popt = scipy_linregress(postshift[time_col]-postshift[time_col].min(),
                                   np.log(postshift[od_col]))

    # Extract the stall biomass and the postshift initial biomass
    mean_stall_biomass = data[data['phase']=='stall'][od_col].mean()


    # Compute the heuristic lag time
    t_regrowth = postshift['time_hr'].values[0] + (np.log(mean_stall_biomass) - postshift_popt[1])/postshift_popt[0]
    t_stall = data[data['phase']=='stall'][time_col].values[0]
    t_lag = t_regrowth - t_stall

    _df = pd.DataFrame({'preshift_growth_rate':preshift_popt[0],
                        'postshift_growth_rate':postshift_popt[0],
                        'entry_time':t_stall,
                        'exit_time':t_regrowth,
                        'lag_time':t_lag},
                        index=[0])
    return _df



def compute_acclimation_time(species_df, 
                             nutrient_df, 
                             quantity='tRNA_c',
                             nutrient_label=0):
    """
    Compute the acclimation time for charged tRNA in a given ecosystem. This
    corresponds to the maximum derivative of the charged tRNA concentration.

    Parameters
    ----------
    species_df : DataFrame
        A DataFrame containing species data. It must have columns 'time_hr' and 'M'.
    nutrient_df : DataFrame
        A DataFrame containing nutrient data. It must have columns 'nutrient_label', 'env_conc', and 'time_hr'.
    nutrient_label : int
        The label of the nutrient to consider.

    Returns
    -------
    float
        The acclimation time.

    """
    # find the time at which the desired nutrient is exhausted
    t_exhaust = nutrient_df[(nutrient_df['nutrient_label']==nutrient_label) & 
                            (nutrient_df['env_conc'] <= 0)]['time_hr'].values[0]
    M_stall = species_df[species_df['time_hr']==t_exhaust]['M'].values[0]

    # Compute the acclimation time.
    postshift = species_df[species_df['M'] >= M_stall]
    ind_acc = np.argmax(np.diff(postshift[quantity].values))
    t_acc = postshift['time_hr'].values[ind_acc+1] - postshift['time_hr'].values[0]
    return t_acc


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

    # Iterate through each species and dilution phase
    for g, d in species_df.groupby('species_label'):
        n_ss_total = 0
        for _g, _d in d.groupby('dilution_cycle'):
            lam = _d['gamma'] * _d['ribosome_content']
            frac_dlam = np.abs(np.diff(lam)/lam[:-1])
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
