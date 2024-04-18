#%%
import numpy as np 
import pandas as pd 

# Set the experimental constants
DATE = '2024-04-17'
MAPPER = {'M0':{'strain':'WT', 'replicate':1,
                'preshift_carbon':'glucose', 'preshift_conc_M':0.005,
                'postshift_carbon': 'acetate', 'postshift_conc_M':0.030},
          'M1':{'strain':'WT', 'replicate':2,
                'preshift_carbon':'glucose', 'preshift_conc_M':0.005,
                'postshift_carbon': 'acetate', 'postshift_conc_M':0.03},
          'M2':{'strain':'WT', 'replicate':3,
                'preshift_carbon':'glucose', 'preshift_conc_M':0.005,
                'postshift_carbon': 'acetate', 'postshift_conc_M':0.03}}


# Iterate through each device
growth_curves = pd.DataFrame([])
for device, details in MAPPER.items():
    data = pd.read_csv(f'./raw/{DATE}_{device}_main_data.csv')

    # Restrict the keys and rename
    data = data[['exp_time', 'od_measured', 'media_temp',
                 'stirring_rate']]
    data['exp_time'] *= 1/3600
    data.rename(columns={'exp_time':'time_hr',
                         'od_measured':'od650nm',
                         'media_temp':'culture_temp_C',
                         'stirring_rate':'stir_rate'},
                         inplace=True)
    # Add identifiers
    data['date'] = DATE
    for k, v in details.items():
        data[k] = v
    growth_curves = pd.concat([growth_curves, data], sort=False)


# Save the processed growth curves to the file.
growth_curves.to_csv(f'./processed/{DATE}_growth_curves.csv', index=False)