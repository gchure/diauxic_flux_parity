#%%

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.viz
import diaux.model
import diaux.quant
import importlib
importlib.reload(diaux.model)
cor, pal = diaux.viz.matplotlib_style()

# Define the species which can eat two nutrients
suballocation = {'strategy': 'dynamic',
                 'nu_max': [10.0, 2.0],
                 'frac_useful': [1, 1],
                 'K': [1E-5, 1E-5],
                 'n': [1, 1]}
species = diaux.model.FluxParityAllocator(suballocation, label=1)                   

# With the species in place, set up the ecosystem with defined nutrient parameters
# and preculture in the pre-shift condition
nutrients = {'init_concs': [0.001, 4]}#, 0.03]}
ecosystem = diaux.model.Ecosystem([species], nutrients)
ecosystem.preculture(equil_time=100)

#%%
# Grow the ecosystem
species_df, nutrient_df = ecosystem.grow(20, dt=1/60)

#%%
# Find the times where nutrient becomes exhausted
exhaustion_times = pd.DataFrame()
for g, d in nutrient_df.groupby(['nutrient_label', 'dilution_cycle']):
    ind = np.where(d['env_conc']<=0)[0]
    if len(ind) > 0:
        t_exhaust = d.iloc[ind[0]]['time_hr']
        _df = pd.DataFrame({'nutrient_label': g[0],
                        'dilution_cycle': g[1],
                        'exhaustion_time': t_exhaust},
                        index = [0])
        exhaustion_times = pd.concat([exhaustion_times, _df], sort=False)

#%%

#%%
# Compute the exhaustion yield


            # val = ss_v[idx]
# ss_df['duration'] = ss_df['t_stop'] - ss_df['t_start']

#%%
def draft_profile_steady_states(species, species_df):
    M_INIT = species_df['M'].values[0]
    M_STAR= (nutrients['init_concs'][0] * species.Y[0] + M_INIT)
    t_ind = np.where(species_df['M'] >= M_STAR[0])[0][0]
    t_hit = species_df.iloc[t_ind]['time_hr']

    # Identify the steady-state regimes
    alloc_stability_thresh = 5E-3
    tRNA_stability_thresh = 1E-3
    alloc_ss = np.abs(1 - species_df['alloc_stability']) <= alloc_stability_thresh
    tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= tRNA_stability_thresh
    species_df['steady_state'] = alloc_ss * tRNA_c_ss
    ss_df = pd.DataFrame([])
    for g, d in species_df.groupby('species_label'):
        ss_v = d['steady_state'].values
        _diff = np.diff(ss_v)
        inds = np.where(_diff != False)[0]
        if len(inds) > 0:
            start, stop = [], []
            for i, idx in enumerate(inds):  
                if ss_v[idx+1] == False:
                    if (i==0):
                      start.append(0)
                    stop.append(idx)
                else:
                    start.append(idx)
            if len(start) > len(stop):
                stop.append(len(ss_v))

            for i, (_start, _stop) in enumerate(zip(start, stop)):
                __start = species_df.iloc[_start]
                __stop = species_df.iloc[_stop]
                phase = d.iloc[_start:_stop]
                phase_lam = phase['gamma'].values * phase['ribosome_content'].values
                _df = pd.DataFrame({'time_start':__start['time_hr'],
                                    'time_stop':__stop['time_hr'],
                                    'M_start':__start['M'],
                                    'M_stop':__stop['M'],
                                    'growth_rate_hr': phase_lam.mean(),
                                    'species_label': g,
                                    'steady_state_idx': i},
                                    index=[0]) 
                ss_df = pd.concat([_df, ss_df], sort=False)
    return ss_df 
#%%
# Find the lag time
postshift = ss_df[ss_df['steady_state_idx']==1]
lag_time = ((np.log(M_STAR) - np.log(postshift['M_start']))/postshift['growth_rate_hr']) + postshift['time_start']
#%%
plt.plot(species_df['time_hr'], species_df['M'])
plt.hlines(M_STAR, 0, 4, color=cor['primary_red'], lw=1)
plt.vlines(t_hit, 1E15, 1E18, color=cor['primary_blue'], lw=1)
plt.fill_betweenx(np.logspace(15, 18), t_hit, t_hit+lag_time, alpha=0.25)
plt.xlim([0, 15])

plt.yscale('log')
#%% Profile the steady-states
alloc_stability_thresh = 5E-3
tRNA_stability_thresh = 1E-3
alloc_ss = np.abs(1 - species_df['alloc_stability']) <= alloc_stability_thresh
tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= tRNA_stability_thresh
species_df['steady_state'] = alloc_ss * tRNA_c_ss
for g, d in species_df.groupby('species_label'):
    _diff = np.diff(d['steady_state'])
    inds = np.where(_diff != False)

#%%

#%%
# Find the steady-states
alloc_ss = np.abs(1 - species_df['alloc_stability']) <= alloc_stability_thresh
tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= tRNA_stability_thresh
# tRNA_u_ss = np.abs(1 - species_df['tRNA_u_stability']) <= 0.003
lam = species_df['gamma'] * species_df['ribosome_content']
dlam = np.abs(np.diff(lam))
dlam = np.insert(dlam, 0, np.inf)

species_df['steady_state'] = tRNA_c_ss  * (dlam <= 0.0001)
plt.plot(species_df['time_hr'], species_df['steady_state'])
# plt.ylim([-0.01, 0.01])
# plt.yscale('log')

#%%
# for g, d in species_df.groupby(['species', 'dilution_cycle']):
#     nuts = nutrient_df[nutrient_df['dilution_cycle']==1]

#%%
plt.plot(species_df['time_hr'], species_df['M']/ diaux.model.OD_CONV)
plt.yscale('log')
#%%
plt.plot(species_df['time_hr'], species_df['tRNA_c'])


#%%
# Find the time and mass of balance growth
lam = species_df['gamma'] * species_df['ribosome_content']
asp_lam = species_df['gamma'] * species_df['phi_Rb']
plt.plot(species_df['time_hr'], lam, 'b-')
plt.plot(species_df['time_hr'], asp_lam, 'k-')
# plt.plot(species_df['time_hr'], species_df['phi_Rb'] * species_df['gamma'])
# plt.plot(species_df['time_hr'], species_df['phi_Rb'])
#%% 

# Find the steady-states
a_thresh = 5E-3
t_thresh = 1E-3 
alloc_ss = np.abs(1 - species_df['alloc_stability'])
tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability'])
tRNA_u_ss = np.abs(1 - species_df['tRNA_u_stability'])
fig, ax = plt.subplots(3, 1, figsize=(6, 4), sharex=True)
for a in ax:
    a.set_yscale('log')
    a.set_ylim([1E-7, 1.1])
ax[0].plot(species_df['time_hr'], alloc_ss)
ax[1].plot(species_df['time_hr'], tRNA_u_ss)

ax[0].hlines(a_thresh, 0, 8, color=cor['primary_red'])
ax[1].hlines(t_thresh, 0, 8, color=cor['primary_red'])
# ax[2].plot(species_df['time_hr'], species_df['gamma'])
# ax[2].hlines(a_thresh, 0, 8, color=cor['primary_red'])
# 
#%%
a_thresh = 5E-3
t_thresh = 1E-3
alloc_ss = np.abs(1 - species_df['alloc_stability']) <= a_thresh
tRNA_c_ss = np.abs(1 - species_df['tRNA_c_stability']) <= t_thresh
species_df['steady_state'] = alloc_ss * tRNA_c_ss
plt.plot(species_df['time_hr'], species_df['steady_state'])
#%%
M_INIT = 0.04 * diaux.model.OD_CONV
M_STAR= (nutrients['init_concs'][0] * species.Y[0] + M_INIT)
lam_1 = species_df.iloc[3]['gamma'] * species_df.iloc[3]['phi_Rb']
t_hit = np.log(M_STAR/M_INIT)/lam_1
t_hit
#%%

plt.hlines(M_STAR/diaux.model.OD_CONV, 0, 20, color=cor['primary_red'],
           linewidth=2)
plt.vlines(t_hit, 0, 1, color=cor['primary_blue'], lw=2)
plt.plot(species_df['time_hr'], species_df['M'] / diaux.model.OD_CONV, lw=2)
plt.yscale('log')
# plt.yscale('log')
