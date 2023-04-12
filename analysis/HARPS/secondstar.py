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
from scipy import stats
from astropy.io import fits as pyfits
from matplotlib.collections import LineCollection
import scipy.optimize as opt

eso = Eso()
eso.login("RishabhDave")
table = eso.query_surveys('HARPS', cache=False, target="HIP87937")
table.pprint()
dp_id = 'ADP.2017-07-25T01:01:03.749'
data_files = eso.retrieve_data(dp_id)

hdulist = pyfits.open('ADP.2017-07-25T01 01 03.749.fits')
''''''''''
scidata = hdulist[1].data

wave = scidata[0][0]
flux = scidata[0][1]
arr2 = scidata[0][2]

plt.plot(wave, flux)



def find_psf(df, column_name):
    """
    This function takes a Pandas DataFrame and the name of a column containing
    the data to analyze, generates a histogram, and finds the point spread
    function (PSF) based on the full width at half maximum (FWHM).
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The DataFrame containing the data to analyze
    column_name : str
        The name of the column containing the data to analyze
    
    Returns:
    --------
    float
        The point spread function (PSF) based on the full width at half maximum (FWHM)
    """
    
    # Extract the data to analyze
    data = df[column_name].values
    
    # Generate the histogram of the data
    hist, bin_edges = np.histogram(data, bins=100)
    
    # Find the bin with the maximum value
    max_bin_index = np.argmax(hist)
    
    # Find the bin with the half maximum value
    half_max = hist[max_bin_index] / 2.0
    half_max_bins = np.where(hist >= half_max)[0]
    
    # Find the first and last bins with the half maximum value
    first_half_max_bin = half_max_bins[0]
    last_half_max_bin = half_max_bins[-1]
    
    # Calculate the full width at half maximum (FWHM) of the histogram
    fwhm = bin_edges[last_half_max_bin] - bin_edges[first_half_max_bin]
    
    # Calculate the point spread function (PSF) as half the FWHM
    psf = fwhm / 2.0
    
    # Plot the histogram
    plt.hist(data, bins=100, alpha=0.5)
    plt.axvline(x=bin_edges[first_half_max_bin], color='r')
    plt.axvline(x=bin_edges[last_half_max_bin], color='r')
    plt.title('Histogram of ' + column_name)
    plt.xlabel(column_name)
    plt.ylabel('Count')
    plt.show()
    
    return psf

find_psf(hyades, 'dist_50') '''''''''