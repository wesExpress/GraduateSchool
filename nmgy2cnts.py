import numpy as np
from astropy.io import fits

image_file = input("Enter .fits image: ")
band = input("Enter SDSS filter: ")
outfile = input("Enter output file: ")
print(image_file, band)

hdu_list = fits.open(image_file)
image_data = hdu_list[0].data
hdu_list.close()

filters = ['u','g','r','i','z']

b_vals = [1.4,0.9,1.2,1.8,7.4]
f_vals = [0.0130280,0.00379653,0.00520193,0.00653877,0.0332798]

ind = filters.index(band)
print(ind)
b = b_vals[ind]
f = f_vals[ind]

cnts = image_data/f

final_image = cnts
hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile, clobber = True)


