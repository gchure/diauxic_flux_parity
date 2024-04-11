#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import diaux.viz
cor, pal = diaux.viz.matplotlib_style()

data = pd.read_csv('../../../data/bioreactor_growth/2024-04-09_pilot_glucose_acetate_shift_dilutions.csv')
data = data[data['od_measured'] > 0]
# data['exp_time'] *= 1/60
plt.plot(data['exp_time'][::5], data['od_measured'][::5], 'k-', lw=1)
# plt.yscale('log')
# plt.xlim([0, 20])
# plt.ylim([1E-3, 0.8])
# plt.yscale('log')