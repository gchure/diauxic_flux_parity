import numpy as np 
import pandas as pd 
from scipy.ndimage import median_filter as scipy_median_filter
from scipy.stats import pearsonr as scipy_pearsonr

def annotate_diauxic_phases(time, 
                            optical_density,
                            pearson_thresh=0.97,
                            pearson_window=30,
                            median_filter=True,
                            median_filter_window=11):
    """
    Classifies portions of a diauxic shift growth curve using a Pearson-correlation
    rolling filter. 

    Parameters
    ----------
    time : numpy.ndarray 
        The time dimension of the growth curve. 
    optical_density : numpy.ndarray
        The  measured optical density of the the growth curve as a function of time.
    pearson_thresh : float (0, 1]
        The threshold above which the log optical density is considered to be 
        linear and, thus, in the exponential growth phase. 
    pearson_window: integer
        The size of the window over which to compute the correlation coefficient. 
        Default is 30 points
    median_filter : bool
        If True, the growth curve will be median filtered before the correlation 
        coefficient is applied. Default is True 
    median_filter_window : odd integer
        The size of the kernel for the median filtering. Only applied if
        `median_filter` is True. Default is 11. 
    
    Returns
    -------
    annotation : pandas DataFrame
        A pandas DataFrame with trimmed data. Regions are labeled as `preshift-exponential`,
        `lag`, and `postshift-exponential`. 
    """
    # Log transform the the optical density and median filt if desired. 
    log_od = np.log(optical_density)
    if median_filter:
        log_od = scipy_median_filter(log_od, size=median_filter_window,
                                             mode='nearest')
        log_label = '_filtered'
    else:
        log_label = ''

    # Compute the pearson correlation coefficients and apply the threshold     
    corr = np.empty(len(log_od) - pearson_window)
    for i in range(len(log_od)-pearson_window):
        corr[i] = scipy_pearsonr(time[i:i+pearson_window], log_od[i:i+pearson_window])[0]
    corr_thresh = corr >= pearson_thresh
    
    # Clip the time and log OD given the window size to avoid indexing slips.
    _time = time[:-pearson_window]
    _log_od = log_od[:-pearson_window]
    _optical_density = optical_density[:-pearson_window]

    # Identify the indices of the transition and apply logic to make sure 
    # the entry/exit logic is in phase. 
    transitions = np.where(np.diff(np.sign(corr_thresh.astype(int))) != 0)[0]
    if len(transitions) == 0:
        raise ValueError("No distinct phases identified!")
    if corr_thresh[0] == False: 
        _time = _time[transitions[0]:] 
        _optical_density = _optical_density[transitions[0]:]
        _log_od = _log_od[transitions[0]:]
        corr = corr[transitions[0]:]

        # Modify the transitions in place to ignore the first non-exponential
        # measurements
        transitions -= transitions[0]
        transitions = transitions[1:]
    if corr_thresh[-1] == False:
        _time = _time[:transitions[-1]]
        _optical_density = _optical_density[:transitions[-1]] 
        corr = corr[:transitions[1]]

    # Create the phase labels
    _phases = np.ones_like(_time)
    _phases[0:transitions[0]] = 0 
    _phases[transitions[1]:] = 2
    mapper = {0:'preshift-exponential',
              1:'lag',
              2:'postshift-exponential'}
    phases = [mapper[p] for p in _phases]
    
    # Pack into a pandas DataFrame and assemble. 
    annotation = pd.DataFrame(np.array([_time, 
                                        _optical_density, 
                                        _log_od,
                                        corr]).T, 
                               columns=['time_hr', 
                                        'optical_density',
                                        'log_optical_density' + log_label,
                                        'correlation_coeff'])
    annotation['annotation'] = phases
                                        
    return annotation

