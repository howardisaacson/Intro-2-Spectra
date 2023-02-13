import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astroquery.eso import Eso
import astropy.wcs as fitswcs #wcs
from specutils import Spectrum1D, SpectralRegion #spectrum1D (specutils)
from pyvo.dal import tap
from astropy import units as u #units
import scipy.interpolate
from scipy.interpolate import splev, splrep
import sys
from astropy.io import fits as pyfits
from matplotlib.collections import LineCollection

eso = Eso()

dp_id = 'ADP.2021-05-30T01:02:35.349'
table = eso.query_surveys('HARPS', cache=False, target="HD10700")
#

hdulist = pyfits.open( "ADP.2021-05-30T01 02 35.343.fits" )

# print column information
hdulist[1].columns

# get to the data part (in extension 1)
scidata = hdulist[1].data

wave = scidata[0][0]
flux = scidata[0][1]
arr2 = scidata[0][2]
# etc.
# where arr1 will contain the data corresponding to the column named: hdulist[1].columns[1]
# where arr2 will contain the data corresponding to the column named: hdulist[1].columns[2]
# etc.

# To plot using maptplotlib:




plt.set_cmap('Spectral')

bin = 100

# this list will contain the indices corresponding to each of the 95th percentile flux values in each bin
indices = []

for i in np.arange((len(wave) - (bin)), step = bin):
    flux_values = []
    for j in np.arange(i, i + bin, step = 1):
        value = flux[j]
        flux_values = np.append(flux_values, value)
    # find the 95th percentile flux value: we use 95 to get the maximum flux value in general 
    # but avoid issues with cosmic rays and other emission lines
    flux_in_bin = np.percentile(flux_values, 95)
    # find the closest value in the flux array to the 95th percentile value
    absolute_difference_function = lambda list_value : abs(list_value - flux_in_bin)
    flux_in_bin = min(flux_values.tolist(), key=absolute_difference_function)
    index_in_bin = flux_values.tolist().index(flux_in_bin)
    index = i + index_in_bin
    indices = np.append(indices, index)

# these lists will contain the wavlength and flux values at each index in 'indices'
wave_values = []
fl_values = []

for index in indices:
    ind = int(index)
    wave_values = np.append(wave_values, wave[ind])
    fl_values = np.append(fl_values, flux[ind])

spl = splrep(wave_values, fl_values, s = 500000)
flux_fit = splev(wave, spl)


first_normalized_flux = flux / flux_fit

plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([5885, 5900])
plt.ylabel('Flux')
plt.title('Na')
plt.show()


plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([4380, 4385])
plt.ylabel('Flux')
plt.title('Fe 1')
plt.show()

plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([4300, 4315])
plt.ylabel('Flux')
plt.title('Ca/Fe')
plt.show()

plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([4850, 4870])
plt.ylabel('Flux')
plt.title('H-β')
plt.show()

plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([5260, 5280])
plt.ylabel('Flux')
plt.title('Fe 2')
plt.show()

plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([6560, 6570])
plt.ylabel('Flux')
plt.title('H-\u03B1')
plt.show()

plt.plot(wave, first_normalized_flux)
plt.xlabel('Wavelength [A]')
plt.xlim([6860, 6870])
plt.ylabel('Flux')
plt.title('Telluric O₂')
plt.show()