#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import tqdm
import diaux.model
import diaux.quant 
import diaux.viz
cor, pal = diaux.viz.matplotlib_style()


# Define the suballocation details for the species
nu_max_preshift = 2
nu_max_postshift = 1
step = 0.1
useful_frac = np.arange(step, 1 + step, step) 
suballocation = {'strategy': 'hierarchical',
                 'K': [1E-4, 1E-4],
                 'n': [1, 1],
                 'Y': [3E19, 1E19],
                 'hierarchy': [0, 1],
                 'nu_max': [nu_max_preshift, nu_max_postshift]}
nutrients = {'init_concs': [0.001, 1.0]}
#%%
species_df = pd.DataFrame()
lag_times = pd.DataFrame()
for i, f in enumerate(tqdm.tqdm(useful_frac)):
    suballocation['frac_useful'] = [1, f]
    species = diaux.model.FluxParityAllocator(suballocation, metabolic_hierarchy=True)
    ecosystem = diaux.model.Ecosystem([species], nutrients)
    ecosystem.preculture(equil_time=100, verbose=False)
    _species_df, _nutrient_df = ecosystem.grow(150, dt=1/60, verbose=False)
    steadystates = diaux.quant.profile_steady_states(_species_df)
    lag, nexh = diaux.quant.draft_profile_lag_time(_nutrient_df, _species_df, steadystates)
    _species_df['lam'] = _species_df['gamma'] * _species_df['ribosome_content']
    tRNAc_tacc = diaux.quant.compute_acclimation_time(_species_df, _nutrient_df, steadystates)
    lam_tacc = diaux.quant.compute_acclimation_time(_species_df, _nutrient_df, quantity='lam')
    _species_df['nu_max_preshift'] = nu_max_preshift
    _species_df['nu_max_postshift'] = nu_max_postshift
    _species_df['useful_frac'] = f
    post_state = _species_df[(_species_df['time_hr'] >= steadystates['time_start'].values[1]) & (_species_df['time_hr'] <= steadystates['time_stop'].values[1])].mean()
    _lag_times = pd.DataFrame({'lag_time_hr': lag['lag_time_hr'].values[0],
                               'preshift_growth_rate': steadystates['avg_growth_rate_hr'].values[0],
                               'postshift_growth_rate': steadystates['avg_growth_rate_hr'].values[1],
                               'tRNA_acclimation_time': tRNAc_tacc,
                               'lam_acclimation_time': lam_tacc,
                               'nu_max_preshift': nu_max_preshift,
                               'nu_max_postshift': nu_max_postshift,
                               'phi_x_pre':  0,
                               'phi_x_post': (1 - f) * post_state['phi_Mb_2'],
                               'useful_frac': f}, index=[0])
    lag_times = pd.concat([lag_times, _lag_times], sort=False)
    species_df = pd.concat([species_df, _species_df], sort=False)

#%% 
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_ylim([0, 15])
ax.set_xlim([-0.05, 0.35])
ax.plot(lag_times['phi_x_post'], lag_times['lag_time_hr'], '-o')
ax.set_xlabel('total useless expression $\phi_{X}$', fontsize=6)
ax.set_ylabel('lag time [hr]', fontsize=6)
# plt.savefig('../../figures/presentations/relative_useless_expression.pdf')
# 
#%%
import seaborn as sns
fig, ax = plt.subplots(1,1, figsize=(3, 2))
ax.set_ylim([0.1, 1])
ax.set_xlim([0, 20])
mapper = {k:v for k, v in zip(useful_frac, sns.color_palette('flare', n_colors=len(useful_frac)))}
for g, d in species_df.groupby('useful_frac'):
    ax.plot(d['time_hr'].values[::5], d['M'].values[::5]/diaux.model.OD_CONV, '-',
            color=mapper[g], lw=1)
ax.set_yscale('log')
ax.set_xlabel('time [hr]', fontsize=6)
ax.set_ylabel('approximate OD$_{600nm}$', fontsize=6)
plt.savefig('../../figures/presentations/relative_useless_expression_curves.pdf')


#%%
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_ylim([0, 15])
ax.set_xlim([-0.05, 0.35])
for g, d in lag_times.groupby('useful_frac'):
    ax.plot(d['phi_x_post'], d['lag_time_hr'], 'o',
            color=mapper[g])

ax.set_xlabel('total useless expression $\phi_{X}$', fontsize=6)
ax.set_ylabel('lag time [hr]', fontsize=6)
plt.savefig('../../figures/presentations/relative_useless_expression.pdf')

#%%