'''
Created on Jul 1, 2015

@author: Sushant S. Mahajan
'''
# Script to read the header of a fits file

from astropy.io import fits

data,header=fits.getdata("1.fits",header=True)
#print header

header['OBJECT']='SUN'
header['TYPE-OBS']='FULLDISK'
header['SEEING']=0
header['CDELT1']=
header['CDELT2']=
header['CTYPE1']='ARCSEC'
header['CTYPE2']='ARCSEC'
header['CRVAL1']=
header['CRVAL2']=
header['TIME-OBS']=
header['DATE-OBS']=
header['DATE_OBS']=
header['ORIGIN']='BBSO'
header['TELESCOPE']='SINGER'
header['WAVELNTH']='HALPHA'
header['OBSERVER']=
header['EXPTIME']=
header['TEC_TEMP']=
# COMMENT ORIGIN FILENAME 
header['CRPIX1']=
header['CRPIX2']=
header['ASP']=
header['CENX']=
header['CENY']=
header['MAXC']=
header['WIDT']=

#fits.writeto('1_bbsohead.fits',data,header,clobber=True)

print header
