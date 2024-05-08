#%%
import numpy as np 
import pandas as pd 
import futileprot.viz
import futileprot.fitderiv
import altair as alt
import tqdm
import altair_saver
import scipy.stats
import scipy.signal
colors, palette = futileprot.viz.altair_style()

# Add metadata
DATE = '2021-09-14'
RUN_NO = 1
STRAINS = 'DoubleKO'
MEDIUM = 'glucose-acetate'
med1, med2 = MEDIUM.split('-')

# Load the measurement data
data = pd.read_csv(f'./output/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_labeled_regions.csv')

info = {k:[[], [], []] for k in data['strain'].unique()}
info_df = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby(['strain', 'replicate'])):

        # Separate by phases 
        med1_phase = d[d['phase']==f'exponential_{med1}'] 
        med2_phase = d[d['phase']==f'exponential_{med2}']

        # Compute the median OD of the lag phase
        lag_od = np.log(d[d['phase']=='shift']['od_600nm_subtracted'].median())
                
        # Infer the parameters for the growth phases using a simple linear regression
        med1_popt = scipy.stats.linregress(med1_phase['elapsed_time_hr'].values,
                                            np.log(med1_phase['od_600nm_subtracted'].values))
        med2_popt = scipy.stats.linregress(med2_phase['elapsed_time_hr'].values,
                                            np.log(med2_phase['od_600nm_subtracted'].values))

        # Compute the two time points
        shift_begin = (lag_od -  med1_popt[1]) / med1_popt[0]
        shift_end = (lag_od - med2_popt[1]) / med2_popt[0]

        # Compute the lag time
        lag_time = shift_end - shift_begin

        # Assemble the dataframe
        info_df = info_df.append({
                        'strain': g[0],
                        'replicate': int(g[1]),
                        f'{med1}_growth_rate': med1_popt[0],
                        f'{med2}_growth_rate':med2_popt[0],
                        'shift_begin': shift_begin,
                        'shift_end': shift_end,
                        'lag_time': lag_time,
                        'lag_od': np.exp(lag_od)
        }, ignore_index=True)

# %%
# Compute the average of the stats
avg_df =  info_df.groupby(['strain']).agg(('mean', 'std')).reset_index()
avg_df
# %%
plot = False
for g, d in data.groupby(['strain']):
        row = False
        for _g, _d in d.groupby(['replicate']):
                _info_df = info_df[(info_df['strain']==g) & 
                                               (info_df['replicate']==_g)]
                t_lag = _info_df.lag_time.values[0]
                data_points = alt.Chart(_d).mark_point().encode(
                        x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                        y=alt.Y('od_600nm_subtracted:Q', title='optical density [a.u.]',
                                scale=alt.Scale(type='log')),
                        color=alt.Color('phase:N')
                )
                info_base = alt.Chart(_info_df)
                shift_begin = info_base.mark_rule().encode(
                                x=alt.X('shift_begin:Q'))
                shift_end= info_base.mark_rule().encode(
                                x=alt.X('shift_end:Q')) 
                lag_od = info_base.mark_rule().encode(
                                y=alt.Y('lag_od:Q')
                        
                ).properties(title=f'{g} rep. {_g} : δ = {t_lag:0.2f} hr')
                layout = data_points + shift_begin + shift_end + lag_od
                if row == False:
                        row = layout
                else:          
                        row |= layout
        if plot == False:
                plot = row
        else:
                plot &= row

altair_saver.save(plot, f'./output/{DATE}_r{RUN_NO}_{STRAINS}_{MEDIUM}_shift_plot.png')

# %%
