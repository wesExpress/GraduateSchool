from astropy.io import fits
import numpy as np
import scipy as sp
import sys as sys

images_file = sys.argv[1]
text_file = sys.argv[2]

hdu_list = fits.open(images_file)
image_data = hdu_list[0].data
hdu_list.close()

np.savetxt(text_file, image_data)
