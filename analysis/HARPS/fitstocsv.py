from astropy.io import fits as pyfits
import csv
from matplotlib import pyplot as plt


# Load the FITS file
hdulist = pyfits.open('ADP.2017-07-25T01 01 03.749.fits')
scidata = hdulist[1].data 

wave = scidata[0][0]
flux = scidata[0][1]
arr2 = scidata[0][2]

plt.plot(wave, flux)

from astropy.io import fits
import csv

# Load the FITS file
fits_file = fits.open('ADP.2017-07-25T01 01 03.749.fits')
data = fits_file[1].data  # Assuming the data is in the first extension

# Get the header keywords and values
header = fits_file[1].header
keywords = header.keys()
values = [header[k] for k in keywords]

# Create the output CSV file
output_file = 'output2.csv'
with open(output_file, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # Write the header keywords and values to the CSV file
    writer.writerow(['Keyword', 'Value'])
    for i, k in enumerate(keywords):
        writer.writerow([k, values[i]])

    # Write the data to the CSV file
    writer.writerow([])  # Add a blank row for readability
    writer.writerow(list(data.columns.names))
    for row in data:
        writer.writerow(list(row))

# Close the FITS file
fits_file.close()
#In this modified version, we use the list() function to convert each row of data to a
#  list of values, which is then written to the CSV file using writer.writerow(). This en
# sures that each data point is written in a separate column below the headers.





