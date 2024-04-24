#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import diaux.viz
cor, pal = diaux.viz.matplotlib_style()

data = pd.read_csv('../../data/bioreactor/pilot_diauxic_shifts.csv')
