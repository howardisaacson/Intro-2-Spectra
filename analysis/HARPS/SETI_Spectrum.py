import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import axes
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
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import random as rd

#eso = Eso()

#dp_id = 'ADP.2021-05-30T01:02:35.349'
#table = eso.query_surveys('HARPS', cache=False, target="HD10700")
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

#plt.scatter(wave, flux, s=1)



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
weights = []


for index in indices:
    ind = int(index)
    wave_values = np.append(wave_values, wave[ind])
    fl_values = np.append(fl_values, flux[ind])
    


spl = splrep(wave_values, fl_values, s = 5000000)
flux_fit = splev(wave, spl)


first_normalized_flux = flux / flux_fit
flux98 = np.percentile(first_normalized_flux, 98)
normalized_flux = first_normalized_flux / flux98
normalized_flux -= 1 + 30 + 289
#plt.plot(wave, normalized_flux)

peak_array = []

def emission_periodic_gaussian(x, amp, mean, stdev):
    '''a is amplitude, b is  the mean (where i think it is), c is standard deviation (sqrt(variance))''' 
    #offset = x[0] #good estimate
    return amp*np.e**(-(x-mean)**2/(2*stdev**2)) + 1



#print(type(fit_maybe))
peak_count = 0
for mean in range(int(wave[0]), int(wave[-1]), 100):
    for _ in range(10):
        peak = emission_periodic_gaussian(wave, rd.randint(0, 200) / 1000, mean, 2.6)
        normalized_flux += peak
        peak_array.append(normalized_flux)
        peak_count += 1
    plt.plot(wave, peak)
    plt.xlim(wave[0], wave[-1])

#print(type(peak))

wave2 = wave
flux2 = normalized_flux
                                              

'''''def gaussian(x, amplitude, mean, sigma):
        #sigma should be the width of point-spread function (PSF)
        #return amplitude * np.exp(-0.5 * ((x - mean) / sigma)**2) + 1
        plt.plot(amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2)))

def laser_inj(wave, flux, wave_spacing, std, amp):
    
    mean_arr = np.arange(wave[0], wave[-1], wave_spacing)
    flux_ = np.copy(flux)
    for i in range(len(mean_arr)):
        flux_ +=  gaussian(wave, amp, mean_arr[i], std) -1
        
    return mean_arr, flux_ '''''''''


plt.plot(wave, normalized_flux)
plt.title('Gaussian Injection')
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux')

axs= plt
plt.show()

def find_peaks_above_threshold(wave2, flux2, threshold):
    """Identify all peaks above a certain threshold in a pyplot object."""
    # Extract x and y data from the plot
    x = wave2
    y = flux2

    # Use SciPy's find_peaks function to identify peaks above the threshold
    peaks, _ = find_peaks(y, height=threshold, distance=2.6)

    # Plot the identified peaks as red dots
    plt.plot(wave[peaks], y[peaks], ms=20, color='red')
    plt.ylim(0, 2.5)

    # Return the identified peaks as a NumPy array
    return peaks

# Example usage

peaks = find_peaks_above_threshold(wave, normalized_flux, 1.1)
#print(type(peaks))

maxima_count = 0

def find_maxima(x_vals, y_vals, threshold, min_distance):
    """
    Finds all local maxima on a graph that exceed a certain threshold value.
    
    Args:
    - x_vals: a list or array of x values
    - y_vals: a list or array of y values
    - threshold: a float representing the threshold value above which a maximum is considered
    
    Returns:
    - maxima: a list of tuples, where each tuple contains the x and y values of a maximum
    """
    maxima = []
    last_max_pos = -min_distance
    for i in range(1, len(y_vals)-1):
        if y_vals[i] > y_vals[i-1] and y_vals[i] > y_vals[i+1] and y_vals[i] > threshold and x_vals[i] - last_max_pos >= min_distance:
            maxima.append((x_vals[i], y_vals[i]))
            last_max_pos = x_vals[i]
    return maxima

# Find all local maxima above threshold value of 0.5
maxima = find_maxima(wave, normalized_flux, 1.1, 1)

def find_minima(wave,flux, amp_threshold,wave_threshold):
    min_arr = np.array([], dtype = int)
    
    under_threshold = np.where(flux < amp_threshold)[0]
    j = 3
    
    for i in under_threshold:
        try:
            if np.all(flux[i - j: i + j+1] >= flux[i]):
                if len(min_arr) == 0:
                    min_arr = np.append(min_arr,i)
                if wave[i] - wave[min_arr[-1]] > wave_threshold:
                    min_arr = np.append(min_arr, i)
        except:
            continue
        
    return min_arr


def laser_rec(wavelength, flux, amp_threshold, wave_threshold):
    
    maxima = find_minima(wavelength,-flux, -amp_threshold,wave_threshold)
    return  wavelength[maxima],flux[maxima]

wave_, flux_ = laser_rec(wave, normalized_flux, 1.1, 4)
plt.plot(wave, normalized_flux)
plt.plot(wave_, flux_, '.', color ='r', label ='recovery')
alteff = len(flux_)
plt.show()

# Plot the data and maxima
hist_flux = []
plt.plot(wave, normalized_flux)
for max in maxima:
    plt.scatter(max[0], max[1], color='red')
    hist_flux.append(max[1])
    maxima_count += 1
plt.show()

''''''''''''''''

local_max = argrelextrema(normalized_flux, np.greater)
plt.plot(wave, normalized_flux)
for i in range(len(normalized_flux)):
    value = normalized_flux[i]
    #print(value)
    if value in local_max[0]:
        plt.scatter(wave[i], normalized_flux[i])
plt.show()

'''''''''''

efficiency =100 * maxima_count / peak_count 
print(f"The efficiency rate of recovery is {efficiency}%")
efficiency2 =100 * alteff / peak_count 
print(f"The alternative efficiency rate of recovery is {efficiency2}%")
# create a histogram of the data
n, bins, patches = plt.hist(hist_flux, bins=4, density=False, alpha=0.75)

# calculate the median value of the data
median = np.median(normalized_flux)

# add a vertical line to the plot at the median value
#plt.axvline(median, color='r', linestyle='dashed', linewidth=1)

# add labels and title to the plot
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of Data')

# show the plot
plt.show()

# print the median value
print(f"The median value is {median}")

plt.show()
recovery = plt
 




#gaussian(flux, 1, 6000, 2.6)
#plt.xlim(5700, 6300)

'''''''''''

def gaussian_curve(x, width, amplitude):
    """Generate a Gaussian curve with given width and amplitude."""
    return amplitude * np.exp(-(x / width) ** 2)


def repeating_gaussian_curve(x, num_repeats, width, amplitude, shift):
    """Generate a repeating Gaussian curve with given number of repeats, width, amplitude, and horizontal shift."""
    x = np.asarray(x)
    y = gaussian_curve(x - shift, width, amplitude)
    for i in range(1, num_repeats):
        y += gaussian_curve(x - shift - i * len(x), width, amplitude)
    return y

'''''''''''
''''''''''

plt.plot(wave, flux, label = 'Data')
plt.scatter(wave_values, fl_values, color = 'black', label = 'Flux Values in the 95th Percentile')
plt.title('Mapping out the Echelle Blaze Function Fit')
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux')
plt.legend()
plt.show()



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

'''''''''''
''''def find_peaks_above_threshold(plot, threshold):
    """Identify all peaks above a certain threshold in a pyplot object."""
    # Extract x and y data from the plot
    lines = axes.Axes.lines(plot)
    x = lines[0].get_xdata()
    y = lines[0].get_ydata()

    # Use SciPy's find_peaks function to identify peaks above the threshold
    peaks, _ = find_peaks(y, height=threshold)

    # Plot the identified peaks as red dots
    plt.plot(x[peaks], y[peaks], 'ro')

    # Return the identified peaks as a NumPy array
    return peaks'''''

# Example usage

y = normalized_flux
plt.plot(wave, y)
peaks = find_peaks_above_threshold(wave, y, threshold=1.1)
plt.show()

print("Identified peaks:", peaks)

#histograms of peaks recovered and injected
#fwhm
