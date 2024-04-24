import numpy as np 
import pandas as pd 
from scipy.signal import savgol_filter as scipy_savgol
from scipy.stats import pearsonr as scipy_pearsonr

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
                          exponential_pearson_thresh=0.98,
                          stall_pearson_thresh=0.75,
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

    # Unpack the correlation to a vector to avoid pandas slowness
    corr = data[corr_col].values

    # Apply the thresholds
    stall = corr <= stall_pearson_thresh
    exp = corr >= exponential_pearson_thresh   
    exp_switchpoints = np.where(np.diff(np.sign(exp.astype(float))) != 0)[0]
    stall_switchpoints = np.where(np.diff(np.sign(stall.astype(float))) != 0)[0]


    # Pad the exponential switchpoints
    exp_switchpoints = np.insert(exp_switchpoints, 0, 0)
    exp_switchpoints = np.append(exp_switchpoints, len(corr)-1)     

    # Compute the phase lengths
    lengths = np.diff(exp_switchpoints)

    # Determine if the data begins in the exponential phase
    exp_phase = exp[0]

    _x = [0]
    for i, ell in enumerate(lengths):
        if ell > minimum_exponential_length:
            _x.append(exp_switchpoints[i+1])
        
    print(exp_switchpoints, lengths, _x)







