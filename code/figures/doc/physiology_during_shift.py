#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.quant
import diaux.model
import diaux.viz 
cor, pal = diaux.viz.matplotlib_style()

# Set the simulation parameters
suballocation = {'strategy':'hierarchical',
                 'nu_max': [3.4, 1.4],
                 'K': [1E-5, 1E-5],
                 'n': [1, 1],
                 'Y': [3E19, 1E19],
                 'hierarchy': [0, 1]}
nutrients = {'init_concs':[0.0005, 0.030]}
species = diaux.model.FluxParityAllocator(suballocation)
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(od_init=0.001, equil_time=100)
species_df, nutrient_df = ecosystem.grow(20, dt=1/60)
steadystates = diaux.quant.profile_steady_states(species_df)
lagtimes, nexh = diaux.quant.draft_profile_lag_time(nutrient_df, species_df, 
                                                    steadystates)
lagtimes
#%%
steadystates
#%%
# Set up a plot to characterize all of the physiology and environmental 
# details
fig, ax = plt.subplots(3, 3, figsize=(6, 4))
ax = ax.ravel()
for a in ax:
    a.set_xlabel('time [hr]', fontsize=6)
ax[8].set_ylim([0, 8])

# ##############
#  FORMATTING  #
# ##############
# Add appropriate labels
ax[0].set_ylabel('$c_{nt_i}(t)/c_{nt_i, 0}$\nrelative nutrient conc.', 
                 fontsize=6)
ax[0].set_title('(A) environment nutrients', fontsize=6)
ax[1].set_title('(B) ecosystem biomass', fontsize=6)
ax[1].set_ylabel('$M(t)/M_0$\nrelative biomass', fontsize=6)
ax[1].set_yscale('log')

ax[2].set_title('(C) instantaneous growth rate', fontsize=6)
ax[2].set_ylabel('$\lambda$ [hr$^{-1}$]\n growth rate',
                fontsize=6)

ax[3].set_title('(D) ribosomes', fontsize=6)
ax[4].set_title('(E) total metabolism', fontsize=6)
ax[5].set_title('(F) metabolic pathway 1', fontsize=6)
ax[6].set_title('(G) metabolic pathway 2', fontsize=6)
for i in range(3, 7):
    ax[i].set_ylabel('resource allocation', fontsize=6)

ax[7].set_title('(H) tRNA pools', fontsize=6)
ax[7].set_ylabel('$tRNA_x/K_{M_x}$\nrelative tRNA concentration', fontsize=6)
ax[8].set_title('(I) translation rate', fontsize=6)
ax[8].set_ylabel('$\gamma(t)$ [hr$^{-1}$]\ntranslation rate', fontsize=6)

# ##############
#  POPULATION  #
# ##############
# Plot the nutrients
for i in range(2):
    if i == 0:
        c = cor['red']
    else:
        c = cor['light_red'] 
    label = f'nutrient {i+1}'
    n = nutrient_df[nutrient_df['nutrient_label']==i]
    ax[0].plot(n['time_hr'], n['env_conc']/n['env_conc'].values[0], '-',
               lw=1, color=c, label=label)

# Plot the species
ax[1].plot(species_df['time_hr'], species_df['M']/species_df['M'].values[0],
           '-', color=cor['primary_black'], lw=1)


# Plot the growth rate
species_df['inst_lam'] = species_df['gamma'] * species_df['ribosome_content']
ax[2].plot(species_df['time_hr'], species_df['inst_lam']/species_df['inst_lam'].values[0],
           '-', lw=1, color=cor['primary_blue'])

# Plot the ribosomal allocation info.
ax[3].plot(species_df['time_hr'], species_df['phi_Rb'], '-', color=cor['primary_gold'],
           lw=1, label='$\phi_{Rb}$')
ax[3].plot(species_df['time_hr'], species_df['ribosome_content'], '--', color=cor['gold'],
           lw=1, label='$M_{Rb}/M$')

# Plot the total metabolic allocation
ax[4].plot(species_df['time_hr'], species_df['phi_Mb_1'] + species_df['phi_Mb_2'], '-', color=cor['primary_purple'],
           lw=1, label='$\phi_{Mb_{tot}}$')
ax[4].plot(species_df['time_hr'], (species_df['M_Mb_1'] + species_df['M_Mb_2'])/species_df['M'], '--', color=cor['purple'],
           lw=1, label='$M_{Mb_{tot}}/M$')

# Plot the metabolic pathway 1 allocation
ax[5].plot(species_df['time_hr'], species_df['phi_Mb_1'], '-', color=cor['primary_purple'],
           lw=1, label='$\phi_{Mb_1}$')
ax[5].plot(species_df['time_hr'], species_df['M_Mb_1']/species_df['M'], '--', color=cor['purple'],
           lw=1, label='$M_{Mb_1}/M$')

# Plot the metabolic pathway 2 allocation
ax[6].plot(species_df['time_hr'], species_df['phi_Mb_2'], '-', color=cor['primary_purple'],
           lw=1, label='$\phi_{Mb_2}$')
ax[6].plot(species_df['time_hr'], species_df['M_Mb_2']/species_df['M'], '--', color=cor['purple'],
           lw=1, label='$M_{Mb_2}/M$')

# Plot the charged and uncharged tRNA dynamics
ax[7].plot(species_df['time_hr'], species_df['tRNA_c']/species.Km_c, '-', color=cor['primary_green'],
           lw=1, label='tRNA$_c$')
ax[7].plot(species_df['time_hr'], species_df['tRNA_u']/species.Km_u, '--', color=cor['green'],
           lw=1, label='tRNA$_u$')

# Plot the translation rates
ax[8].plot(species_df['time_hr'], species_df['gamma'], '-', color=cor['primary_red'],
           lw=1)
# Add legends, where necessary
ax[0].legend(fontsize=5, handlelength=0.5, loc='center')
ax[3].legend(fontsize=5)#, handlelength=1)
ax[4].legend(fontsize=5)#, handlelength=1)
ax[5].legend(fontsize=5)#, handlelength=1)
ax[6].legend(fontsize=5)#, handlelength=1)
ax[7].legend(fontsize=5, handlelength=0.5, loc='center right')

plt.tight_layout()
plt.savefig('../../figures/doc/physiological_shift_dynamics.pdf', bbox_inches='tight')