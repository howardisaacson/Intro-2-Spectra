import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from astroquery.eso import Eso
import astropy.wcs as fitswcs #wcs
from specutils import Spectrum1D, SpectralRegion #spectrum1D (specutils)
from pyvo.dal import tap
from astropy import units as u #units
import astropy.wcs as fitswcs #wcs
from specutils import Spectrum1D, SpectralRegion #spectrum1D (specutils)
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
arr1 = scidata[0][1]
arr2 = scidata[0][2]
# etc.
# where arr1 will contain the data corresponding to the column named: hdulist[1].columns[1]
# where arr2 will contain the data corresponding to the column named: hdulist[1].columns[2]
# etc.

# To plot using maptplotlib:



plt.plot(wave, arr1)

plt.set_cmap('Spectral')

plt.title("HD10700 Spectrum")

plt.xlabel("Wave")

plt.ylabel("Flux")

x    = wave
y    = arr1
cols = np.linspace(0,1,len(x))

plt.show()

