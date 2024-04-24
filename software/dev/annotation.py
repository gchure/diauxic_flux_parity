#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import importlib 
import diaux.data
import diaux.viz
import scipy.signal
import scipy.stats
importlib.reload(diaux.data)
cor, pal = diaux.viz.matplotlib_style()

data1 = pd.read_csv('./data/2024-04-22_complete_growth_curves.csv')
data2 = pd.read_csv('./data/2024-04-17_complete_growth_curves.csv')
data2 = data2[data2['tech_rep'] != 3]
valid_reps = [7, 2, 3]
data1 = data1[data1['tech_rep'].isin(valid_reps)]
data =  pd.concat([data1, data2], sort=False)
data['idx'] = data.groupby(['date', 'tech_rep']).ngroup()
data = data[(data['od650nm'] >= 0.1) &  (data['od650nm'] <= 1.1)]

corr_dfs = []
for g, d in data.groupby('idx'):
    _corr_df = diaux.data.compute_pearson_correlation(d)
    corr_dfs.append(_corr_df)
corr_df = pd.concat(corr_dfs, sort=False)

_corr_df = corr_df[corr_df['idx'] == 3]
diaux.data.classify_diauxic_phases(_corr_df)

#%%
# Step 1: savgol filter
dfs = []
for g, d in data.groupby('idx'):
    d = d.copy()
    d['filtered'] = scipy.signal.savgol_filter(np.log(d['od650nm']),
                                               window_length=30,
                                               polyorder=1)
    dfs.append(d)
filt_data = pd.concat(dfs, sort=False)

fig, ax = plt.subplots(2, 3, figsize=(6, 4), sharex=True)
ax = ax.ravel()
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel('log optical density', fontsize=6)
for g, d in data.groupby('idx'):
    ax[g].plot(d['time_hr'], np.log(d['od650nm']), '.', ms=3, markeredgewidth=0)
    ax[g].set_title(f'replicate {g}', fontsize=6)

for g, d in filt_data.groupby('idx'):
    ax[g].plot(d['time_hr'], d['filtered'], 'r-', lw=1)
ax[-1].axis(False)


#%%
# Compute the pearson correlation coefficient with a length 30 window. 
dfs = []
pearson_window = 30
for g, d in filt_data.groupby('idx'):
    corr = np.zeros(len(d) - pearson_window)
    for i in range(len(corr)):
        corr[i] = scipy.stats.pearsonr(d['time_hr'].values[i:i+pearson_window], 
                                       d['filtered'].values[i:i+pearson_window])[0]
    d = d[:-pearson_window].copy()
    d['pearson_r'] = corr
    dfs.append(d)
corr_df = pd.concat(dfs, sort=False)

fig, ax = plt.subplots(2, 3, figsize=(6,4), sharex=True)
ax = ax.ravel()
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel('Pearson correlation coeff.', fontsize=6)
for g, d in corr_df.groupby('idx'):
    zero_crossings = np.where(np.diff(np.sign(d['pearson_r'].values)) != 0)[0]
    # ax[g].fill_betweenx([-1, 1], d['time_hr'].values[zero_crossings[0]],
                        # d['time_hr'].values[zero_crossings[1]], 
                        # alpha=0.5, color='k')
    ax[g].plot(d['time_hr'], d['pearson_r'], 'k-')
    ax[g].plot(d['time_hr'], np.sign(d['pearson_r']), 'r-')
ax[-1].axis(False)
# for a in ax:
    # a.set_yscale('log')

#%%
# Label exponential regions (corr > exp_thresh) and stalled regions based on 
# zero crossings
exp_thresh = 0.995
stall_thresh = 0.5
dfs = []
for g, d in corr_df.groupby('idx'):
    d = d.copy()
    zero_crossings = np.where(np.diff(np.sign(d['pearson_r'].values)) != 0)[0]
    t_start = d['time_hr'].values[zero_crossings[0]-1]
    t_stop = d['time_hr'].values[zero_crossings[1]-1]
    d['label'] = 'transition' 
    d.loc[(d['pearson_r'] >= exp_thresh) & (d['time_hr'].values < t_start),
          'label'] = 'preshift exponential'
    d.loc[(d['pearson_r'] >= exp_thresh) & (d['time_hr'].values > t_stop),
          'label'] = 'postshift exponential'
    d.loc[(d['time_hr'].values >= t_start) & 
          (d['time_hr'].values <= t_stop), 'label'] = 'stall'
    dfs.append(d)
label_df = pd.concat(dfs, sort=False)

fig, ax = plt.subplots(2, 3, figsize=(6, 4))
ax = ax.ravel()
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel('Pearson correlation coeff.', fontsize=6)

for g, d in label_df.groupby(['idx', 'label']):
    ax[g[0]].plot(d['time_hr'], np.log(d['od650nm']), '.', ms=3, markeredgewidth=0,
                  label=g[1])

ax[-1].axis(False)
ax[0].legend()

#%%
filt_data = pd.concat(dfs, sort=False)

fig, ax = plt.subplots(2, 3, figsize=(6, 4), sharex=True)
ax = ax.ravel()
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel('inst. growth rate [hr$^{-1}$]', fontsize=6)

for g, d in filt_data.groupby('idx'):
    ax[g].plot(d['time_hr'].values[1:], np.diff(d['filtered'].values)/np.diff(d['time_hr'].values), 'r-', lw=1)
ax[-1].axis(False)

plt.tight_layout()