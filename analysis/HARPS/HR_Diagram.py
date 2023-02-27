import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()
import astropy.units as u
import astropy.coordinates as coord
## making a GAIA cone_search of 20m radius around M45 center

#Now we can read the csv file into a pandas dataframe

hyades = pd.read_csv('hyades.csv') # I renamed my csv file to 'm5.csv' and put it in the the subfolder 'data'
hyades['bp_rp'] = hyades['phot_bp_mean_mag'] - hyades['phot_rp_mean_mag']
hyades['abs_mag'] = hyades['phot_g_mean_mag'] - 5 * np.log10((hyades['dist_50']*0.01)) + 5

t = list(hyades['dist_50'])
sns.scatterplot(x = 'bp_rp', y='phot_g_mean_mag', data = hyades, c = t[0:26998], cmap = 'RdBu', s = 5)

plt.ylabel('Absolute Magnitude (G band)')
plt.xlabel('B-R Color')
plt.title('Hyades Cluster HR Diagram')
#plt.xlim(0,3)
#sns.scatterplot(hr.b_v+0.5, hr.V+10)

plt.gca().invert_yaxis()
plt.show()

polaris = pd.read_csv('polaris.csv') # I renamed my csv file to 'm5.csv' and put it in the the subfolder 'data'
polaris['bp_rp'] = polaris['phot_bp_mean_mag'] - polaris['phot_rp_mean_mag']
polaris['abs_mag'] = polaris['phot_g_mean_mag'] - 5 * np.log10((polaris['dist_50']*0.01)) + 5

t = list(polaris['dist_50'])
sns.scatterplot(x = 'bp_rp', y='phot_g_mean_mag', data = polaris, c = t[0:21152], cmap = 'RdBu', s = 5)

plt.ylabel('Absolute Magnitude (G band)')
plt.xlabel('B-R Color')
plt.title('Polaris HR Diagram')
#plt.xlim(0,3)
#sns.scatterplot(hr.b_v+0.5, hr.V+10)

plt.gca().invert_yaxis()
plt.show()

M3 = pd.read_csv('M3.csv') # I renamed my csv file to 'm5.csv' and put it in the the subfolder 'data'
M3['bp_rp'] = M3['phot_bp_mean_mag'] - M3['phot_rp_mean_mag']
M3['abs_mag'] = M3['phot_g_mean_mag'] - 5 * np.log10((M3['dist_50']*0.01)) + 5

t = list(M3['dist_50'])
sns.scatterplot(x = 'bp_rp', y='phot_g_mean_mag', data = M3, c = t[0:2411], cmap = 'RdBu', s = 5)

plt.ylabel('Absolute Magnitude (G band)')
plt.xlabel('B-R Color')
plt.title('M3 Cluster HR Diagram')
#plt.xlim(0,3)
#sns.scatterplot(hr.b_v+0.5, hr.V+10)

plt.gca().invert_yaxis()
plt.show()

galacticcenter = pd.read_csv('galacticcenter.csv') # I renamed my csv file to 'm5.csv' and put it in the the subfolder 'data'
galacticcenter['bp_rp'] = galacticcenter['phot_bp_mean_mag'] - galacticcenter['phot_rp_mean_mag']
galacticcenter['abs_mag'] = galacticcenter['phot_g_mean_mag'] - 5 * np.log10((galacticcenter['dist_50']*0.01)) + 5

t = list(galacticcenter['dist_50'])
sns.scatterplot(x = 'bp_rp', y='phot_g_mean_mag', data = galacticcenter, c = t[0:2760], cmap = 'RdBu', s = 5)

plt.ylabel('Absolute Magnitude (G band)')
plt.xlabel('B-R Color')
plt.title('Galactic Center HR Diagram')
#plt.xlim(0,3)
#sns.scatterplot(hr.b_v+0.5, hr.V+10)

plt.gca().invert_yaxis()
plt.show()