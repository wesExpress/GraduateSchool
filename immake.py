from astropy.io import fits
from astropy.table import Table
import numpy as np
import sys as sys

print ''
print 'Reading in text file...'
im_file = sys.argv[1]
im_fits = sys.argv[2]
im_img = np.loadtxt(im_file)
print 'Done.'
print 'Writing image...'
fits.writeto(im_fits, im_img)
print 'Done.'
print ''

