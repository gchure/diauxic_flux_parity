#%%
import numpy as np 
import pandas as pd 
import diaux.model
import diaux.quant

# Define experiment parameters
DATE = '2021-09-08'
STRAINS = 'SingleKO'
MEDIUM = 'glucose-acetate'
RUN_NO = 1
ROOT = '../../../..'
SKIPROWS = 36 

# Identify the mdia
med1, med2 = MEDIUM.split('-')
OD_BOUNDS = [0.01, 0.15]
ROLLING_WINDOW = 8 # Number of points to consider for rolling correlation window
NUDGE = 2
PEARSON_THRESH = 0.975 # Threshold for pearson correlation

# Add the well identifiers
MAP = {'GC032': ['C3', 'D3', 'E3'],
       'GC049': ['C4', 'D4', 'E4'],
       'GC052': ['C5', 'D5', 'E5'],
       'GC047': ['C6', 'D6', 'E6'],
       'GC050': ['C7', 'D7', 'E7'],
       'GC048': ['C8', 'D8', 'E8'],
       'GC053': ['C9', 'D9', 'E9'],
       'GC055': ['C10', 'D10' ,'E10'],
       'GC030': ['F3', 'F4', 'F5'],
       'GC029': ['F6', 'F7', 'F8'],
       'GC001': ['F9', 'F10', 'F11']} 

# Generate a list of all valid wells
wells = [f'{letter}{number}' for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] for number in np.arange(1,13)]

# Load the data
data = pd.read_csv(f'./raw/{DATE}_r{RUN_NO}.csv', 
                skiprows=SKIPROWS)

# Melt and drop unnecessary stuff
melted = data.melt(id_vars=['Time'], var_name='well', value_name='od600nm')
melted = melted.loc[melted['well'].isin(wells)]
melted.dropna(inplace=True)

# Add strain identifier and replicates
melted['strain'] = 'blank'
melted['replicate'] = 0
for strain, wells in MAP.items():
    for idx, well in enumerate(wells):
        melted.loc[melted['well']==well, 'strain'] = strain
        melted.loc[melted['well']==well, 'replicate'] = idx + 1

# Add information regarding date and growth medium
melted['growth_medium'] = MEDIUM
melted['date'] = DATE
melted['run_number'] = RUN_NO

# Convert time to elapsed time
melted['time_sec'] = pd.to_timedelta(melted['Time'].values)
melted['time_sec'] = melted['time_sec'].dt.total_seconds()
melted['time_hr'] = (melted['time_sec'] - melted['time_sec'].min())/3600

# Drop unnecessary Time columns
melted.drop(columns=['Time', 'time_sec'], inplace=True)


# Reformat blank value as average eentry per time
measurement = []
for g, d in melted.groupby(['time_hr']):
    avg_blank = d[d['strain']=='blank']
    meas = d[d['strain']!='blank'].copy()
    meas['avg_blank_value'] = avg_blank['od600nm'].mean()
    measurement.append(meas)
measurement = pd.concat(measurement, sort=False)
measurement.rename(columns={'strain':'identifier'}, inplace=True)

# Add shorthand strain information and class identifier
strain_shorthand, _, strain_class = diaux.io.standardize_strains(measurement['identifier'].values)
measurement['strain'] = strain_shorthand
measurement['class'] = strain_class

# Perform the blank subtraction
measurement['od600nm_subtracted'] = measurement['od600nm'].values - measurement['avg_blank_value'].values


# Save to disk
measurement.to_csv(f'./processed/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_measurements.csv', index=False)


#%%
# Perform the labeling of the regions
trunc = measurement[(measurement['od600nm_subtracted'] >= OD_BOUNDS[0]) & 
                    (measurement['od600nm_subtracted'] <= OD_BOUNDS[1])]

phases = pd.DataFrame([])
lags = pd.DataFrame([])
for g, d in trunc.groupby(['strain', 'replicate', 'class']):
        # Rescale the time to the minimum
        d['time_hr'] -= d['time_hr'].min()
        corr_df = diaux.quant.compute_pearson_correlation(d,
                                                          time_col='time_hr',
                                                          od_col='od600nm_subtracted',
                                                          savgol_filter=False,
                                                          pearson_window=8)
        _phases = diaux.quant.classify_diauxic_phases(corr_df,
                                             od_col='od600nm_subtracted',
                                             exponential_pearson_thresh=PEARSON_THRESH,
                                             minimum_exponential_length=5)
        _phases['medium'] = MEDIUM
        _lag = diaux.quant.quantify_lag_time(_phases, od_col='od600nm_subtracted')
        _lag['date'] = DATE
        _lag['run_number'] = RUN_NO
        _lag['strain'] = g[0]
        _lag['replicate'] = g[1]
        _lag['class'] = g[2]
        _lag['medium'] = MEDIUM
        lags = pd.concat([lags, _lag], sort=False)
        phases = pd.concat([phases, _phases], sort=False)
 

phases.to_csv(f'./processed/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_labeled_phases.csv', index=False)
lags.to_csv(f'./processed/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_lag_times.csv', index=False)