#%%
import pandas as pd 
import glob 
files = glob.glob('./2024*/processed/*valid_shifts.csv')
df = pd.concat([pd.read_csv(f) for f in files], sort=False)
df.to_csv('../../../data/bioreactor/pilot_diauxic_shifts.csv', index=False)