#%%
import numpy as np 
import pandas as pd 
import glob 
import diaux.model
import diaux.quant 
import diaux.viz
cor, pal = diaux.viz.matplotlib_style()

files = glob.glob('../experiments/plate_reader/*SingleKO*/processed/')