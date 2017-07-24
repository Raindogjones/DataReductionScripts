from astropy.io import fits
import os
import glob


hdulist=fits.open("subtest.fit",mode='update')
hdulist[0].header['CDELT1']=10.0
hdulist[0].header['CRVAL1']=0
hdulist[0].header['CRUNIT1']="km/s"
hdulist[0].header['CRPIX1']=1.0
hdulist[0].header['CTYPE1']="VELO-WAV"

# hdulist[0].header['CD2_3']=0.00004815
# hdulist[0].header['CRVAL2']=281.458
# hdulist[0].header['CRPIX2']=21.0
# hdulist[0].header['CTYPE2']="RA---TAN"
# 
# hdulist[0].header['CD3_3']=0.00004815
# hdulist[0].header['CRVAL3']=-33.3422
# hdulist[0].header['CRPIX3']=1.0
# hdulist[0].header['CTYPE3']="DEC--TAN"

#hdulist[0].header['CD2_2']=0.00004815
hdulist[0].header['CDELT2']=0.52
hdulist[0].header['CRVAL2']=0.0
hdulist[0].header['CRPIX2']=20.0
hdulist[0].header['CRUNIT2']="Arcseconds"
#hdulist[0].header['CTYPE3']="DEC--TAN"

hdulist[0].header['OBJECT']="IC4776"

hdulist.verify('fix')
hdulist.flush()
