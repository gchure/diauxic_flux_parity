#%%
import pandas as pd 
import glob 

# Growth curves
files = glob.glob('./2024*/processed/*valid_shifts_labeled.csv')
df = pd.concat([pd.read_csv(f) for f in files], sort=False)
df.to_csv('../../../data/bioreactor/pilot_diauxic_shifts.csv', index=False)

# Lag times
files = glob.glob('./2024*/processed/*lag_time.csv')
df = pd.concat([pd.read_csv(f) for f in files], sort=False)
df.to_csv('../../../data/bioreactor/pilot_lag_times.csv', index=False)
