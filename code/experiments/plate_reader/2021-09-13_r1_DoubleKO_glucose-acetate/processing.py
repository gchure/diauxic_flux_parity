"""
Notes:

* This script drops several replicates or entire strains due to instrumentation
problems or unknown extreme variance in technical replicates. This script 
drops the following samples:

+ Replicate 1 of ∆opp ∆flh -- Extreme variance from other technical replicates
+ Entirety of ∆his ∆rbs -- Exponential phase growth was not met for post-shift
+ Entirety of ∆opp ∆rbs -- Exponential phase growth was not met for post-shift
+ First point of replicate 1 ∆rbs ∆glt -- First point is obviously an outlier and was dropped.
"""
# %%
import numpy as np 
import pandas as pd 
import futileprot.io
import futileprot.fitderiv
import futileprot.viz
import altair as alt 
import altair_saver
colors, palette = futileprot.viz.altair_style()

# Define experiment parameters
DATE = '2021-09-14'
STRAINS = 'DoubleKO'
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

MAP = {
       'GC064': ['C2', 'D2', 'E2'],
       'GC065': ['C3', 'D3', 'E3'],
       'GC066': ['C4', 'D4', 'E4'],
       'GC067': ['C5', 'D5', 'E5'],
       'GC068': ['C6', 'D6', 'E6'],
       'GC069': ['C7', 'D7', 'E7'],
       'GC070': ['C8', 'D8', 'E8'],
       'GC071': ['C9', 'D9', 'E9'],
       'GC072': ['C10', 'D10' ,'E10'],
       'GC073': ['C11', 'D11' ,'E11'],
       'GC074': ['F3', 'F4', 'F5'],
       'GC075': ['F6', 'F7', 'F8'],
       'GC076': ['F9', 'F10', 'F11']} 

# Generate a list of all valid wells
wells = [f'{letter}{number}' for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] for number in np.arange(1,13)]

# Load the data
data = pd.read_csv(f'{ROOT}/data/diauxic_shifts/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}/{DATE}_r{RUN_NO}.csv', 
                skiprows=SKIPROWS)

# Melt and drop unnecessary stuff
melted = data.melt(id_vars=['Time'], var_name='well', value_name='od_600nm')
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
melted['elapsed_time_hr'] = (melted['time_sec'] - melted['time_sec'].min())/3600

# Drop unnecessary Time columns
melted.drop(columns=['Time', 'time_sec'], inplace=True)


# Reformat blank value as average eentry per time
measurement = []
for g, d in melted.groupby(['elapsed_time_hr']):
    d = d.copy()
    avg_blank = d[d['strain']=='blank']
    meas = d[d['strain']!='blank']
    meas['avg_blank_value'] = avg_blank['od_600nm'].mean()
    measurement.append(meas)
measurement = pd.concat(measurement, sort=False)
measurement.rename(columns={'strain':'identifier'}, inplace=True)

# Add shorthand strain information and class identifier
strain_shorthand, _, strain_class = futileprot.io.standardize_strains(measurement['identifier'].values)
measurement['strain'] = strain_shorthand
measurement['class'] = strain_class

# Perform the blank subtraction
measurement['od_600nm_subtracted'] = measurement['od_600nm'].values - measurement['avg_blank_value'].values

# Drop the noted samples
dfs = []
for g, d in measurement.groupby(['strain']):
    if (g== '∆his ∆rbs') | (g == '∆opp ∆rbs'):
        continue
    else:
        if g == '∆dpp ∆flh':
            dfs.append(d[d['replicate'] != 1])
        elif g == '∆rbs ∆glt':
            for _g, _d in d.groupby(['replicate']):
                if _g == 1:
                    __d = _d[(_d['od_600nm_subtracted'] < 0.02) & (_d['elapsed_time_hr'] < 0.17)]
                    __d['elapsed_time_hr'] -= __d['elapsed_time_hr'].min()
                    dfs.append(__d)
                else:
                    dfs.append(_d)
        else:
            dfs.append(d)

# Save to disk
measurement = pd.concat(dfs, sort=False)
measurement.to_csv(f'./output/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_measurements.csv', index=False)


#%%
# Perform the labeling of the regions
trunc = measurement[(measurement['od_600nm_subtracted'] >= OD_BOUNDS[0]) & 
                    (measurement['od_600nm_subtracted'] <= OD_BOUNDS[1])]



data_dfs = []
gp_dfs = []
for g, d in trunc.groupby(['strain', 'replicate', 'class']):
        # Rescale the time to the minimum
        d['elapsed_time_hr'] -= d['elapsed_time_hr'].min()

        # Using the data, do a rolling correlation coefficient.
        d['log_od'] = np.log(d['od_600nm_subtracted'])
        out = d[['elapsed_time_hr', 
                    'log_od']].rolling(ROLLING_WINDOW).corr().reset_index()
        out = out[out['level_1']=='elapsed_time_hr']['log_od']
        out = out[ROLLING_WINDOW:]
        locs = out >= PEARSON_THRESH
        locs = locs.astype(int)
        sign_loc = np.sign(locs).diff()
        min_ind = np.argmin(sign_loc) + ROLLING_WINDOW - NUDGE
        max_ind = np.argmax(sign_loc) + ROLLING_WINDOW + NUDGE 

        # Convert the min and max ind into times
        exp_med1_time = d['elapsed_time_hr'].values[min_ind]
        exp_med2_time = d['elapsed_time_hr'].values[max_ind]
                        
        # Label the gp and data frames
        d['phase'] = 'shift'
        d.loc[d['elapsed_time_hr'] <=exp_med1_time, 'phase'] = f'exponential_{med1}'
        d.loc[d['elapsed_time_hr'] >=exp_med2_time, 'phase'] = f'exponential_{med2}'

        # append the labeled data frames
        data_dfs.append(d)


data_labeled = pd.concat(data_dfs, sort=False)

data_labeled.to_csv(f'./output/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_labeled_regions.csv', index=False)
#%%


# Generate a figure of all of the raw traces
raw_traces = alt.Chart(
                    data=data_labeled, 
                    width=400, 
                    height=200
                ).mark_point(
                    opacity=0.75
                ).encode(
                    x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                    y=alt.Y('od_600nm_subtracted:Q', title='optical density [a.u.]',
                            scale=alt.Scale(type='log')),
                    color=alt.Color('phase:N', title='identified phase'),
                    shape=alt.Color('replicate:N', title='technical replicate')
                ).facet(
                    row='strain'
                )
altair_saver.save(raw_traces, f'output/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_raw_traces.png',
                 scale_factor=2)
# %%

# %%
