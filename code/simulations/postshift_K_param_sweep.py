#%%
import importlib
import numpy as np 
import pandas as pd
import diaux.model
import diaux.quant
import diaux.viz
import tqdm
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()
K_range = np.logspace(-9, 0, 20)
nu_preshift = 10
nu_postshift = np.linspace(1, nu_preshift, 20)
nutrients = {'init_concs':[0.0005, 30]}
lag_mat = np.zeros((len(K_range), len(nu_postshift)))
lag_df = pd.DataFrame([])
prealloc = np.zeros_like(K_range)
for i, K in enumerate(tqdm.tqdm(K_range)):
    for j, nu_post in enumerate(nu_postshift):
        suballocation = {'strategy': 'hierarchical',
                         'K': [K, 1E-5],
                         'n': [1, 1],
                         'nu_max':[nu_preshift, nu_post]}
        species = diaux.model.FluxParityAllocator(suballocation)
        eco = diaux.model.Ecosystem([species], nutrients)
        eco.preculture(verbose=False, init_conc_override=[1, 0])
        species_df, nutrient_df = eco.grow(50, dt=1E-3, verbose=False)
        steadystates = diaux.quant.profile_steady_states(species_df)
        steadystates['duration'] = steadystates['time_stop'] - steadystates['time_start']
        steadystates = steadystates[steadystates['duration']>=0.1]
        lagtimes, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, steadystates[steadystates['total_ss_idx']==steadystates['total_ss_idx'].max()])       
        if len(lagtimes) > 0:
            if lagtimes['lag_time_hr'].values[0] <= 0.001:
                short = True
            else:
                short = False
        if (len(lagtimes) == 0) | (short==True):
            lagtimes = pd.DataFrame({'dilution_cycle': 1,
                                     'species_label':0, 
                                     'preshift_nutrient_label':0,
                                     'stall_biomass':eco._M0 + species.Y[0] * nutrients['init_concs'][0],
                                     'entry_time_hr':nexh['exhaustion_time_hr'].values[0],
                                     'exit_time_hr':nexh['exhaustion_time_hr'].values[0],
                                     'lag_time_hr':0}, index=[0])

        lagtimes['nu_preshift'] = nu_preshift
        lagtimes['nu_postshift'] = nu_post
        lagtimes['K'] = K
        lagtimes['lam_postshift'] = steadystates['avg_growth_rate_hr'].values[-1]
        lagtimes['lam_preshift'] = steadystates['avg_growth_rate_hr'].values[0]
        lagtimes['alpha_2_preallocation'] = species_df[(species_df['time_hr']>= steadystates.iloc[0]['time_start']) & 
                                                       (species_df['time_hr']<=steadystates.iloc[0]['time_stop'])]['alpha_2'].mean()
        prealloc[i] = lagtimes['alpha_2_preallocation'].values[0]
        lag_df = pd.concat([lag_df, lagtimes], sort=False)
        lag_mat[i, j] = lagtimes['lag_time_hr'].values[0]

#%%
lag_df.to_csv('../../data/simulations/postshift_K_param_sweep_lagtimes.csv')

#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1, subplot_kw={'projection':'3d'}, figsize=(3,3))
X, Y = np.meshgrid(np.log10(K_range), nu_postshift)
# P, Y = np.meshgrid(np.log10(prealloc),)
# _lag_mat = lag_mat + (np.ones_like(lag_mat) * 1E-5) * (lag_mat == 0)
ax.plot_surface(X, Y, np.log10(lag_mat), cmap='viridis')
ax.view_init(20, 30)
ax.set_xticks([0, -2, -4, -6, -8])
ax.set_xticklabels(['10$^0$', '10$^{-2}$', '10$^{-4}$', '10$^{-6}$', '10$^{-8}$'])
ax.set_xlabel('allocation affinity constant\n K$_{A}$ [M]',
              fontsize=6)
ax.set_ylabel('postshift metabolic rate\n' + r'$\nu_{max_2}$', fontsize=6)
ax.set_zlabel('$t_{lag}$ [hr]\nlag time', fontsize=6)
ax.set_zticks([-2, -1, -0, 1])
ax.set_zticklabels(['10$^{-2}$', '10$^{-1}$', '10$^{0}$', '10$^{1}$'])
# ax.set_zlim([-5, 1])

plt.tight_layout()
# ax.set_ylabel(r'')
# import matplotlib.pyplot as plt 
# plt.imshow(lag_mat, origin='lower')

#%%
fig, ax = plt.subplots(1,1)
im = ax.imshow(lag_mat, cmap='viridis', origin='lower')
ax.grid(False)
ax.set_xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
ax.set_xticklabels(np.round(nu_postshift[::2], decimals=1))
ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
ax.set_yticklabels(np.round(np.log10(K_range)[::2], decimals=1))
# cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im,ax=ax, label='lag time [hr]')
# plt.colorbar()