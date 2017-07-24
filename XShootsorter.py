import glob
from astropy.io import fits
import os
for file in glob.iglob('*.fits'):
    hdu=fits.open(file)
    header=hdu[0].header
    if header['HIERARCH ESO PRO CATG']=='SCI_SLIT_FLUX_IDP_VIS':
		os.system("cp %s VIS/." %file)
    if header['HIERARCH ESO PRO CATG']=='SCI_SLIT_FLUX_IDP_UVB':
		os.system("cp %s UVB/." %file)
    if header['HIERARCH ESO PRO CATG']=='SCI_SLIT_FLUX_IDP_NIR':
		os.system("cp %s NIR/." %file)