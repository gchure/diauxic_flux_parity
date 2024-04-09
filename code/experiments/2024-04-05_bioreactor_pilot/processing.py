#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.viz
cor, pal = diaux.viz.matplotlib_style()

data = pd.read_csv('../../../data/bioreactor_growth/2024-04-05_pilot_glucose_acetate_shift.csv')
data = data[data['od_measured'] > 0]
plt.plot(data['exp_time']/2600, data['od_measured'], 'k-', lw=1)
# plt.xlim([0, 20])
# plt.ylim([1E-3, 0.8])
# plt.yscale('log')