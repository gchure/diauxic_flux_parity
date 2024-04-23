#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.viz
cor, pal = diaux.viz.matplotlib_style()

data = pd.read_csv('../../../data/bioreactor_growth/2024-04-10_NCM3722_glucose_dithering.csv')
# data = data[data['od_measured'] > 0]
data['exp_time'] /= 3600 
data = data[(data['od_measured'] >= 0.35) & (data['od_measured'] <= 0.4) & (data['pump_1_rate'] == 0)]
data['exp_time'] -= data['exp_time'].min()

switch_points = np.where(np.round(np.diff(data['od_measured'].values), decimals=2) <=-0.05)[0]
times = data['exp_time'].values
# data.loc[data['exp_time'] <= times[switch_points][1], 'cycle'] = 0
for i, idx in enumerate(switch_points):
    if i == 0:
        pass
    data.loc[(data['exp_time'] > times[switch_points[i-1]+1]) & 
         (data['exp_time'] < times[idx]), 'cycle'] = i
# data.loc[data['exp_time'] >times[idx], 'cycle'] = i+2
# data = data[['exp_time', 'od_measured', 'cycle']] 
# data
# plt.plot(data['exp_time'], data['od_measured'], 'k.', lw=1)
# plt.yscale('log')
# plt.xlim([0, 20])
# plt.ylim([1E-3, 0.8])
# plt.yscale('log')
# np.diff(data['exp_time'])
import scipy.stats
lam = []
for g, d in data.groupby('cycle'):
    
    popt = scipy.stats.linregress(d['exp_time']-d['exp_time'].min(), np.log(d['od_measured']))
    lam.append(popt[0])
    plt.plot(d['exp_time'] -d['exp_time'].min(), d['od_measured'], '-o')

{#%%
n_cols = 3
n_rows = int(np.ceil(data['cycle'].max()/n_cols))
fig, ax = plt.subplots(n_rows, n_cols, figsize=(6, 2 * n_rows), 
                       sharey=True)
_ax = ax.ravel()
for a in _ax:
    a.set_xlabel('time [hr]', fontsize=6)
for g, d in data.groupby('cycle'):
    time = d['exp_time'].values - d['exp_time'].min()
    _ax[int(g-1)].plot(time, np.log(d['od_measured']), 'o', ms=4, color=cor['primary_black'])
    fit = scipy.stats.linregress(time, np.log(d['od_measured']),)
    time_range = np.linspace(0, np.max(time) * 1.02, 10)
 
    _ax[int(g-1)].plot(time_range, fit[1] + fit[0] * time_range, '-',
                       color=cor['primary_red'], lw=1)
    _ax[int(g-1)].set_title(f'Cycle {int(g)}; $\lambda$={fit[0]:0.3f} per hr', 
                            fontsize=6)

_ax[-1].axis('off')
for i in range(n_rows):
    ax[i, 0].set_ylabel('log OD$_{600nm}$', fontsize=6)
plt.tight_layout()
plt.savefig('./2024-04-10_glucose_dithered_growth_fits.pdf', bbox_inches='tight')

#%%
fig, ax = plt.subplots(1,1, figsize=(6, 2))

data = pd.read_csv('../../../data/bioreactor_growth/2024-04-10_NCM3722_glucose_dithering.csv')
# data = data[data['od_measured'] > 0]
data['exp_time'] /= 3600 
data['exp_time'] -= data['exp_time'].min()
ax.plot(data['exp_time'], data['od_measured'], 'k-')
ax.set_ylim([0, 0.61])
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_ylabel('OD$_{600nm}$', fontsize=6)
plt.savefig('./2024-04-19_glucose_dithered_growth_curve.pdf', bbox_inches='tight')
