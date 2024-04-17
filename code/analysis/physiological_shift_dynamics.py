#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.model 
import diaux.viz 
import seaborn as sns
import tqdm
cor, pal = diaux.viz.matplotlib_style()

# Set up storage dataframes for the simulation outputs
species_df = pd.DataFrame([])
nutrient_df = pd.DataFrame([])

# Define postshift ranges
nu_max_preshift = 5
step = 0.5
nu_max_postshift = np.arange(step, nu_max_preshift+step, step)

# Define the common suballocation
suballocation = {'strategy': 'hierarchical',
                 'K': [1E-5, 1E-5],
                 'n': [1, 1]}
nutrients = {'init_concs': [0.0005, 10]}

# Iterate through the shift conditions 
for i, nu_post in enumerate(tqdm.tqdm(nu_max_postshift)):
    # Instantiate and grow the ecosystem
    suballocation['nu_max'] = [nu_max_preshift, nu_post]
    species = diaux.model.FluxParityAllocator(suballocation)
    ecosystem = diaux.model.Ecosystem([species], nutrients)
    ecosystem.preculture(od_init=0.001, equil_time=100, verbose=False)
    _species_df, _nutrient_df = ecosystem.grow(50, dt=1/60, verbose=False)
                                         
    for _d in [_species_df, _nutrient_df]:
        _d['nu_max_preshift'] =  nu_max_preshift
        _d['nu_max_postshift'] = nu_post
    species_df = pd.concat([species_df, _species_df], sort=False)
    nutrient_df = pd.concat([nutrient_df, _nutrient_df], sort=False)


#%%
# Set the color palette
palette = sns.color_palette('crest', n_colors=len(lags_med['nu_max_postshift'].unique()))
cmap = {n:c for n, c in zip(lags_med['nu_max_postshift'].unique(), palette)}

#%%
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
ax = ax.ravel()
ax[0].set_yscale('log')
ax[0].set_xlim([-3, 8])
ax[0].set_ylim([1E-2, 1E3])
ax[1].set_xlim([-3, 8])
ax[1].set_yscale('log')
ax[1].set_ylim([1E-3, 2])
ax[2].set_ylim([0,0.025])
ax[2].set_xlim([0, 8])
for g, d in species_df.groupby('nu_max_postshift'):
    d = d[d['tRNA_c']/d['tRNA_c'].values[0] > 1E-3].copy()
    _nutrients = nutrient_df[nutrient_df['nu_max_postshift']==g]
    nexh  = _nutrients.loc[(_nutrients['nutrient_label']==0) & 
                   (_nutrients['env_conc'] <= 0)]
    t_shift = nexh['time_hr'].values[0]
    d['rescaled_time'] = d['time_hr'] - t_shift
    m_stall = d[d['time_hr']==t_shift]['M'].values[0]
    max_deriv_idx =  np.argmax(np.diff(d['tRNA_c'].values))
    max_deriv_time = d['rescaled_time'].values[max_deriv_idx+1]
    max_deriv_val =  d['tRNA_c'].values[max_deriv_idx+1]/d['tRNA_c'].values[0]
    ax[0].plot(d['time_hr']-t_shift, d['M']/m_stall, '-', color=cmap[g], lw=1)
    ax[1].plot(d['time_hr']-t_shift, d['tRNA_c']/d['tRNA_c'].values[0], color=cmap[g], lw=1)   
    ax[1].plot(max_deriv_time, max_deriv_val, 'o', ms=4, color=cmap[g])
    ax[2].plot(d['time_hr'].values[1:]-t_shift, np.diff(d['tRNA_c'].values)/d['tRNA_c'].values[0], color=cmap[g], lw=1)   
    ax[2].plot(max_deriv_time, np.diff(d['tRNA_c'].values)[max_deriv_idx]/d['tRNA_c'].values[0], 'o',
               ms=4, color=cmap[g])
    ax[3].plot(g, max_deriv_time, 'o', ms=4, color=cmap[g])
#%%
# # Load the simulated data sets
# traj = pd.read_csv('../../data/simulations/preshift_postshift_sweep_species_trajectories.csv')
# lags = pd.read_csv('../../data/simulations/preshift_postshift_sweep_lag_times.csv')

# # Takes care of an edge case where a lag is detected after complete nutrient exhaustion
# lags.dropna(inplace=True) 

# # Select an example case
# traj_med = traj[traj['nu_max_preshift']==5]
# lags_med = lags[lags['nu_max_preshift']==5]

# #%%
# # Set the color palette
# palette = sns.color_palette('crest', n_colors=len(lags_med['nu_max_postshift'].unique()))
# cmap = {n:c for n, c in zip(lags_med['nu_max_postshift'].unique(), palette)}

# #%%
# fig, ax = plt.subplots(2, 2, figsize=(6, 4))
# ax = ax.ravel()
# ax[0].set_yscale('log')
# ax[0].set_xlim([-3, 10])
# ax[0].set_ylim([1E-2, 1E3])
# for g, d in traj_med.groupby('nu_max_postshift'):
#     _lag = lags_med[lags_med['nu_max_postshift']==g]
#     t_shift = _lag['entry_time_hr'].values[0]
#     m_stall = _lag['stall_biomass'].values[0]
#     ax[0].plot(d['time_hr']-t_shift, d['M']/m_stall, '-', color=cmap[g], lw=1)
#     ax[1].plot(d['time_hr']-t_shift, d['tRNA_c']/)   
