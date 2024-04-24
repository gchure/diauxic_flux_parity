import numpy as np 
import pandas as pd 
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