#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model
import diaux.viz 
cor, pal = diaux.viz.matplotlib_style()

# Load the trajectories
traj = pd.read_csv('../../data/simulations/nu_sweep_trajectories.csv')
lag = pd.read_csv('../../data/simulations/nu_sweep_lagtimes.csv')
steadystates = pd.read_csv('../../data/simulations/nu_sweep_steadystates.csv')

#%%
# Select the example 
nu_pre = 10.0
nu_post = 6.0

growth_curve = traj[(traj['nu_max_preshift']==nu_pre) &
                    (traj['nu_max_postshift']==nu_post)]
ss = steadystates[(steadystates['nu_max_preshift']==nu_pre) & 
                  (steadystates['nu_max_postshift']==nu_post)]
_lag = lag[(lag['nu_max_preshift']==nu_pre) & 
           (lag['nu_max_postshift']==nu_post)]
#%%
fig, ax = plt.subplots(1, 1, figsize=(3,2))
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_ylabel('optical density [OD$_{600nm}$/mL]', fontsize=6)
ax.set_yscale('log')
ax.set_xlim([0, 4])
ax.set_ylim([0.04, 10])

ax.plot(growth_curve['time_hr'], growth_curve['M'] / diaux.model.OD_CONV, '-',
        lw=2, color=cor['primary_black'])

# plot the preshift growth rate
preshift_time = np.linspace(0, 8, 200)
preshift_growth = 0.04 * np.exp(ss[ss['steady_state_idx']==0]['growth_rate_hr'].values[0] * preshift_time)
ax.plot(preshift_time, preshift_growth, '--', lw=1, color=cor['primary_red'])
ax.text(1, 0.35, 'pre-shift exponential growth', rotation=42, fontsize=5, color=cor['primary_red'])
# Plot the postshift growth rate
postshift_time = np.linspace(0, 16, 200)
postshift_growth = (_lag['shift_biomass'].values[0]/diaux.model.OD_CONV) * np.exp(postshift_time * ss['growth_rate_hr'].values[0])
ax.plot(postshift_time + _lag['regrowth_time'].values[0], postshift_growth, '--', lw=1, color=cor['primary_green'])
ax.text(2.5, 0.3, 'post-shift exponential growth', rotation=34, fontsize=5, color=cor['primary_green'])

ax.fill_betweenx([0, 10], _lag['shift_time'], _lag['regrowth_time'], color=cor['light_blue'], alpha=0.5, zorder=1)
ax.text(0.95, 0.045, 'lag phase', fontsize=5, color=cor['dark_blue'])

ax.hlines(_lag['shift_biomass'] / diaux.model.OD_CONV, 0, 4, lw=1, color=cor['primary_gold'], linestyle=':')
ax.text(2.7, 0.12, 'biomass at nutrient shift', color=cor['primary_gold'], fontsize=5)
plt.savefig('../../figures/sandbox/example_lag_phase.pdf', bbox_inches='tight')