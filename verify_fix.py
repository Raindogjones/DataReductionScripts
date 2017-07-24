from astropy.io import fits
import glob

for file in glob.iglob('*.fits'):
	hdulist=fits.open(file,mode='update')
	hdulist.verify('fix')
	hdulist.flush()
