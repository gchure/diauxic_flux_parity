#%%
import numpy as np 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import tqdm
import diaux.model
import diaux.viz 
import diaux.quant
import seaborn as sns   
cor, pal = diaux.viz.matplotlib_style()


#%%
nu_max_preshift = 3
nu_max_postshift = np.arange(1, 3.5, 0.1)

suballocation = {'strategy':'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'Y': [3E19, 1E19]}
nutrients = {'init_concs': [0.001, 1.0]}

species_df = pd.DataFrame([])
tacc = pd.DataFrame([])
# Iterate through each nu_max and run the simulation
for i, nu_post in enumerate(tqdm.tqdm(nu_max_postshift)):
    suballocation['nu_max'] = [nu_max_preshift, nu_post]
    species = diaux.model.FluxParityAllocator(suballocation)
    ecosystem = diaux.model.Ecosystem([species], nutrients)
    ecosystem.preculture(equil_time=100, verbose=False)
    _species_df, _nutrient_df = ecosystem.grow(30, dt=1/60, verbose=False)
    steadystates = diaux.quant.profile_steady_states(_species_df)
    _species_df['lam'] = _species_df['gamma'] * _species_df['ribosome_content']
    lagtimes, nexh = diaux.quant.draft_profile_lag_time(_nutrient_df, _species_df, 
                                                        steadystates)
    tRNA_tacc = diaux.quant.compute_acclimation_time(_species_df, _nutrient_df)
    lam_tacc = diaux.quant.compute_acclimation_time(_species_df, _nutrient_df, quantity='lam')

    _tacc = pd.DataFrame({'nu_max_preshift':nu_max_preshift,
                          'nu_max_post': nu_post,
                          'lam_preshift': steadystates['avg_growth_rate_hr'].values[0],
                          'lam_postshift': steadystates['avg_growth_rate_hr'].values[1],
                          'tRNAc_tacc': tRNA_tacc,
                          'lam_tacc': lam_tacc,
                          'lag_time': lagtimes['lag_time_hr'].values[0]}, index=[0])
    _species_df['nu_max_preshift'] = nu_max_preshift
    _species_df['nu_max_postshift'] = nu_post
    species_df = pd.concat([species_df, _species_df], sort=False)
    tacc = pd.concat([tacc, _tacc], sort=False)

#%%
colors = sns.color_palette('crest', n_colors=len(nu_max_postshift))
mapper = {k:v for k, v in zip(nu_max_postshift, colors)}
fig, ax = plt.subplots(1, 3, figsize=(6, 2))
for g, d in species_df.groupby('nu_max_postshift'):
    ax[0].plot(d['time_hr'], d['tRNA_c']/d['tRNA_c'].values[0], color=mapper[g])
    ax[1].plot(d['time_hr'].values[:-1], np.diff(d['tRNA_c'])/d['tRNA_c'].values[0],
               color=mapper[g])
ax[0].set_xlim([1, 10])
ax[1].set_ylim([-0.01, 0.02])
ax[1].set_xlim([2.5, 10])
ax[2].plot(tacc['lag_time'], tacc['tRNAc_tacc'], 'o', color=cor['primary_green'],
         label='tRNA$_c$', ms=4)
ax[2].plot(tacc['lag_time'], tacc['lam_tacc'], 'o', color=cor['primary_blue'],
         label='$\lambda_i$', ms=4)
ax[2].plot([2, 10], [2, 10], '-', color=cor['primary_black'], lw=0.5, label='1-to-1')
ax[2].set_xlabel('lag time [hr]', fontsize=6)
ax[2].set_ylabel('acclimation time [hr]', fontsize=6)
ax[2].legend()
for i in range(2):
    ax[i].set_xlabel('time [hr]', fontsize=6)
ax[0].set_ylabel('relative tRNAc', fontsize=6)
ax[1].set_ylabel('tRNAc derivative', fontsize=6)
ax[2].set_yscale('log')
ax[2].set_xscale('log')
plt.tight_layout()

#%%
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.set_xlabel('lag time [hr]', fontsize=6)
ax.set_ylabel('acclimation time [hr]', fontsize=6)
_id = np.linspace(2, 10, 100)
ax.plot(_id, _id, 'k-', lw=1, color=cor['primary_black'])
ax.plot(tacc['lag_time'], tacc['tRNAc_tacc'], 'o', color=cor['primary_green'], alpha=0.75, 
        label='tRNA acclimation time')
plt.legend()
plt.savefig('../../../figures/presentations/acclimation_time_tRNA.pdf', bbox_inches='tight')
ax.plot(tacc['lag_time'], tacc['lam_tacc'], 'o', color=cor['primary_blue'], alpha=0.75, 
        label='inst. growth acclimation time')
plt.legend()
plt.savefig('../../../figures/presentations/acclimation_time_tRNA_lam.pdf', bbox_inches='tight')


