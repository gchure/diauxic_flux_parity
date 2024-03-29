#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import diaux.viz 
import diaux.model
cor, pal = diaux.viz.matplotlib_style()

# Load an example of the preshift variation
trajectories = pd.read_csv('../../../data/simulations/preshift_variation_trajectories.csv')
steadystates = pd.read_csv('../../../data/simulations/preshift_variation_steadystates.csv')
lagtimes = pd.read_csv('../../../data/simulations/preshift_variation_lagtimes.csv')


#%%
nu_post = 5
nu_pre = 10
_ = [] 
for d in [trajectories, steadystates, lagtimes]:
    _.append(d[(d['nu_max_preshift']==nu_pre) & 
             (d['nu_max_postshift']==nu_post)])
traj, ss, lag = _

#%%
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax[0,0].set_ylabel('$OD_{600nm} / mL$\napproximate optical density', fontsize=6)
ax[0, 1].set_ylabel('$\gamma$ [hr$^{-1}$]\n translation rate', fontsize=6)
ax[1, 0].set_ylabel('$\phi_{Rb}$\n ribosomal allocation', fontsize=6)
ax[1, 1].set_ylabel('$\phi_{Mb}$\n metabolic allocation', fontsize=6)
ax[0,0].set_yscale('log')
ax[0,0].set_ylim([0.04, 5])
for a in ax.ravel():
    a.set_xlabel('time [hr]', fontsize=6)
    a.set_xlim([0, 5])

# Plot the growth curve
ax[0,0].plot(traj['time_hr'], traj['M']/diaux.model.OD_CONV, '-', 
           color=cor['primary_black'], lw=1)

# Plot the translation rate
ax[0,1].plot(traj['time_hr'], traj['gamma']/traj['gamma'].values[0], '-',
           color=cor['primary_blue'], lw=1)
ax[1,0].plot(traj['time_hr'], traj['phi_Rb'], '-',
             color=cor['primary_gold'], lw=1)
ax[1,1].plot(traj['time_hr'], traj['phi_Mb_1'], '-', 
             color=cor['light_purple'], lw=1,
             label='$\phi_{Mb,1}$',zorder=1000)
ax[1,1].plot(traj['time_hr'], traj['phi_Mb_2'], '-', 
             color=cor['dark_purple'], lw=1,
             label='$\phi_{Mb,2}$')
ax[1, 1].legend()
plt.savefig('../../../figures/presentations/diauxic_shift_curves.pdf')

#Add the demarcation of the steady states
preshift = ss[ss['steady_state_idx']==0]
postshift = ss[ss['steady_state_idx']==1]

time_range = np.linspace(0, 5, 200)
ax[0,0].hlines(lag['shift_biomass'].values[0]/diaux.model.OD_CONV, 0, 5, linestyle='--',
               color=cor['light_black'])
ax[0,0].plot(time_range, 0.04 * np.exp(time_range * preshift['growth_rate_hr'].values[0]),
             '--', lw=1, color=cor['primary_red'])
ax[0,0].plot(time_range + lag['regrowth_time'].values[0], 
             (lag['shift_biomass'].values[0]/diaux.model.OD_CONV)*np.exp(time_range * postshift['growth_rate_hr'].values[0]),
             '--', lw=1, color=cor['primary_green'])
ax[0,0].text(0.6, 0.2, 'preshift exponential growth', color=cor['primary_red'],
             fontsize=5, rotation=52)
ax[0,0].text(3, 0.22, 'postshift exponential growth', color=cor['primary_green'],
             fontsize=5, rotation=42)
ax[0,0].text(4, 0.13, 'shift biomass', color=cor['light_black'],
             fontsize=5)
plt.savefig('../../../figures/presentations/diauxic_shift_curves_growth_annotation.pdf')

lims = []
for i, a in enumerate(ax.ravel()):
    lims.append([a.get_ylim()[0], a.get_ylim()[1]])
    a.fill_betweenx(lims[i],
                    lag['shift_time'].values[0], lag['regrowth_time'].values[0],
                    color=cor['pale_blue'], zorder=1, alpha=0.5)
    if i == 0:
        y = 0.05
    else:
        y = 0.5 * lims[i][1]
    a.text(1.25, y, '   $\Delta$ [hr]\nlag time', fontsize=5, color=cor['primary_blue'])
for l, a in zip(lims, ax.ravel()):
    a.set_ylim(l)
plt.savefig('../../../figures/presentations/diauxic_shift_curves_total_annotation.pdf')