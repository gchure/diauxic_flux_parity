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



#%%


corr_dfs = []
for g, d in data.groupby('idx'):
    _corr_df = diaux.data.compute_pearson_correlation(d)
                                                      
    corr_dfs.append(_corr_df)
corr_df = pd.concat(corr_dfs, sort=False)

_corr_df = corr_df[corr_df['idx'] == 4]

#%%

fig, ax = plt.subplots(1,1, figsize=(3, 2))
ax.plot(_corr_df['time_hr'], _corr_df['pearson_correlation'], lw=1)


#%%

lab = diaux.data.classify_diauxic_phases(_corr_df, exponential_pearson_thresh=0.95,
                                         stall_pearson_thresh=0.8)


for g, d in lab.groupby('phase'):
    plt.plot(d['time_hr_centered'], d['od650nm_centered'], 'o', ms=3, markeredgewidth=0,
             label=g)
plt.legend()
# plt.plot(_corr_df['time_hr'], _corr_df['pearson_correlation'])
#%%

lag_df = diaux.data.quantify_lag_time(lab)