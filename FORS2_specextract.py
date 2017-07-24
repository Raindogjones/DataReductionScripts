import glob
from astropy.io import fits
import os
import numpy as np
from PyAstronomy import pyasl


for file in glob.iglob('mapped_sci*.fits'):
    hdu=fits.open(file)
    num=file.replace("mapped_sci_lss_00","")[:-5]
    header=hdu[0].header
    if header['HIERARCH ESO INS SHUT EXPTIME']>900:
    	os.system("/star/bin/convert/./fits2ndf %s wavecal.sdf" %file)
    	ycen=np.argmax(np.sum(hdu[0].data[120:150],axis=1))+120
    	os.system("/star/bin/figaro/./isubset IMAGE=wavecal YSTART=min YEND=max XSTART=min XEND=1950 OUTPUT=trim")
    	print "Extracting spectrum at y=%d" %ycen
    	os.system("/star/bin/figaro/./polysky IMAGE=trim YS1=%d YE1=%d YS2=%d YE2=%d DEGREE=3 NREJECT=2 OUTPUT=sky" %(ycen-13,ycen-8,ycen+8,ycen+13))
    	os.system("/star/bin/figaro/./profile IMAGE=sky YSTART=%d YEND=%d DEGREE=3 NREJECT=2 PROFILE=prof RESIDUAL=res" %(ycen-6,ycen+6))
    	os.system("/star/bin/figaro/./optextract IMAGE=sky PROFILE=prof SPECTRUM=spec")
    	os.system("/star/bin/convert/./ndf2fits spec.sdf spec%s.fits" %num)
    	JD=header["MJD-OBS"]+2400000.5+(0.5*header["HIERARCH ESO INS SHUT EXPTIME"]/86400)
    	helicor,HJD=pyasl.helcorr(header["HIERARCH ESO TEL GEOLON"],header["HIERARCH ESO TEL GEOLAT"],header["HIERARCH ESO TEL GEOELEV"],header["RA"],header["DEC"],JD,debug=False)
    	hdu2=fits.open("spec%s.fits" %num,mode='update')
    	hdu2[0].header["HELICOR"]=helicor
    	hdu2[0].header["HJD"]=HJD
    	hdu2[0].header["EXPTIME"]=header['HIERARCH ESO INS SHUT EXPTIME']
    	hdu2.flush()
    	
os.system("rm *.sdf")
    	
    
