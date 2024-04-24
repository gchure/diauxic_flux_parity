#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.model 
import diaux.viz 
cor, pal = diaux.viz.matplotlib_style()
OD_BOUNDS = [0.1, 1]
VALID_REPS = [2]
DATE = '2024-04-18'
data = pd.read_csv(f'./processed/{DATE}_complete_growth_curves.csv')

# Restrict the data to the valid bounds
data = data[(data['od650nm'] >= OD_BOUNDS[0]) & 
            (data['od650nm'] <= OD_BOUNDS[1])]

# Set up a figure showing all of the replicates
N_COL = 2
N_REP = data['tech_rep'].max()
N_ROWS = int(np.ceil(N_REP / N_COL))
resid =  N_ROWS * N_COL % data['tech_rep'].max()

fig, ax = plt.subplots(N_ROWS, N_COL, figsize=(6, 2 * N_ROWS))
ax = ax.ravel()
for i in range(-resid, 0):
    ax[i].axis(False)
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_ylabel('log optical density', fontsize=6)
for g, d in data.groupby('tech_rep'):
    ax[g-1].plot(d['time_hr'], np.log(d['od650nm']), '.', ms=4, markeredgewidth=0)
    ax[g-1].set_title(f'Device: M{g-1}')

plt.tight_layout()
plt.savefig('./growth_curves.png', bbox_inches='tight', dpi=150)

#%%
data = data[data['tech_rep'].isin(VALID_REPS)]
data.to_csv(f'./processed/{DATE}_valid_shifts.csv', index=False)

data = data[data['tech_rep'].isin(VALID_REPS)]
data.to_csv(f'./processed/{DATE}_valid_shifts.csv', index=False)

labeled_data = pd.DataFrame([])
lag_data = pd.DataFrame([])
for g, d in data.groupby('tech_rep'):
    corr_df = diaux.data.compute_pearson_correlation(d)
    lab_df = diaux.data.classify_diauxic_phases(corr_df,
                                                exponential_pearson_thresh=0.95,
                                                stall_pearson_thresh=0.8)
    lag_df = diaux.data.quantify_lag_time(lab_df)
    lag_df['date'] = DATE
    lag_df['tech_rep'] = g
    lag_df['strain'] = d['strain'].values[0]
    lag_df['preshift_carbon'] = d['preshift_carbon'].values[0]
    lag_df['preshift_conc_M'] = d['preshift_conc_M'].values[0]
    lag_df['postshift_conc_M'] = d['postshift_conc_M'].values[0]
    lag_data = pd.concat([lag_data, lag_df], sort=False)
    labeled_data = pd.concat([labeled_data, lab_df], sort=False)

labeled_data.to_csv(f'./processed/{DATE}_valid_shifts_labeled.csv', index=False)
lag_data.to_csv(f'./processed/{DATE}_lag_time.csv', index=False)

#%%
# Plot the aligned curves
fig, ax = plt.subplots(1,1, figsize=(3, 2))
ax.set_xlabel('time from shift [hr]', fontsize=6)
ax.set_ylabel('relative stall optical density', fontsize=6)
ax.set_yscale('log')
for g, d in labeled_data.groupby('tech_rep'):
    ax.plot(d['time_hr_centered'], d['od650nm_centered'], label=g,
            lw=1)
ax.legend(title='technical replicate')
plt.savefig('./aligned_growth_curves.png', dpi=150)