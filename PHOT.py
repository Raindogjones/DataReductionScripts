from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import alipy
import glob
import Tkinter as tk
import ttk
import os
from os.path import exists
from astropy.time import Time
from time import gmtime, strftime
from astropy import units as u
from astropy.coordinates import SkyCoord
from PyAstronomy import pyasl
from astropy.time import TimeDelta
from math import log10
from math import sqrt

from PythonPhot import getpsf
from PythonPhot import aper
from PythonPhot import pkfit
from PythonPhot import rdpsf

import photutils
from astropy.stats import sigma_clipped_stats
from photutils.psf import create_prf
from photutils.psf import subtract_psf
from photutils.psf import psf_photometry
from photutils.psf import GaussianPSF

import statistics

from photutils import aperture_photometry, CircularAperture, CircularAnnulus
from astropy.table.operations import hstack


#Name one images as ref.fit, how?


images=sorted(glob.glob("EFOSC*.fit"))
images2=sorted(glob.glob("alipy_out/EFOSC*.fits"))
images2=[w.replace('alipy_out/','') for w in images2]
images2=[w.replace('_gregister.fits','.fit') for w in images2]
images_to_align = [item for item in images if item not in images2]
ref_image="ref.fit"

if not exists (ref_image):
	print "No reference image found for alignment and/or star selection"
	quit()

if images_to_align:
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications: 
		if id.ok == True:
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
                
                
	outputshape = alipy.align.shape(ref_image)

	for id in identifications:
		if id.ok == True:
			alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)

LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Helvetica", 10)
SMALL_FONT = ("Helvetica", 8)


def popupmsg(title,msg):
    popup = tk.Tk()
    popup.wm_title(title)
    label = ttk.Label(popup, text=msg, font=LARGE_FONT)
    label.pack(side="top", fill="x", pady=10, padx=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()
    
    
popupmsg("GAIA", "Use Pick Object function to select target and at \nleast two comparisons, saving their coordinates \neach time to GaiaPick.log \nClick okay to open ref.fit in Gaia and begin selection")

# Open ref.fit to select stars

# try:
# 	os.remove('GaiaPick.log')
# except:
# 	pass
# 
# os.system('/star/bin/gaia/gaia.sh ref.fit')

xpos,ypos=np.array([0]),np.array([0])
apcoords=list()

gaialog="GaiaPick.Log"
log=open(gaialog)
for line in log:
	if line.startswith('ref.fit'):
		elements=line.split()
		xpos=np.append(xpos,float(elements[1]))
		ypos=np.append(ypos,float(elements[2]))
		apcoords.append((float(elements[1],float(elements[2]))

xpos=np.delete(xpos,0,0)
ypos=np.delete(ypos,0,0)

coords=np.stack((xpos,ypos),axis=-1)

#stats = open("stats.dat",'w')
os.chdir('alipy_out')
os.system('rm *psf.fits')
differential = open("differential.dat",'w+')
print >>differential, '#HJD \t Tar-c \t Tar-c_err \t comp \t comp_err'


for file in glob.iglob('*.fits'):
	hdu=fits.open(file)
	header=hdu[0].header
	data=hdu[0].data

#Aperture photometry commands - remember to use apcoords (list) and not coords (np.array)	
	mean, median, std =sigma_clipped_stats(data,sigma=5.0)
	threshold = mean + 10.5*std
	sources=photutils.irafstarfind(data,threshold=threshold, fwhm=6.0)
	seeing=statistics.median(sources["fwhm"])
	
	ap_r=1.5*seeing
	an_rin=3*seeing
	an_rout=4.5*seeing
	ap_area=np.pi*ap_r**2
	an_area=np.pi*(an_rin**2 - an_rout **2)
	apertures=CircularAperture(apcoords,r=ap_r)
	annulus_apertures=CircularAnnulus(apcoords, r_in=an_rin, r_out=an_rout)
	rawflux=aperture_photometry(data,apertures)
	bkgflux=aperture_photometry(data, annulus_apertures)
	phot_table=hstack([rawflux,bkgflux],table_names['raw', 'bkg'])
	final_sum=(phot_table['aperture_sum_raw'] - phot_table['aperture_sum_bkg'] * ap_area / an_area)				
	phot_table['residual_aperture_sum'] = final_sum
	starflux=phot_table['residual_aperture_sum'][0:3]
	stardflux=np.sqrt(starflux)
	

#photutils psf photometry commands
# 	mean, median, std =sigma_clipped_stats(data,sigma=5.0)
# 	threshold = mean + 10.5*std
# 	sources=photutils.daofind(data,fwhm=6.0,threshold=threshold)
# 	
# 	stars=np.array([[0,0]])
# 	for source in sources:
# 		stars=np.append(stars,[[source['xcentroid'],source['ycentroid']]],axis=0)
# 
# 	stars=np.delete(stars,[0],axis=0)
# 
# 	prf_discrete = create_prf(data,stars,11,mode='median',subsampling=5)
# 	#fluxes is fluxes of stars with coords from GAIA _not_ from the daophot sources
# 	fluxes=psf_photometry(data,coords,prf_discrete)
# 	starflux=fluxes[0:3]
# 	stardflux=np.sqrt(starflux)



#PythonPhot commands
# 	mag,magerror,flux,fluxerr,sky,skyerr,badflag,outstr = aper.aper(data,xpos,ypos,phpadu=1,apr=5,zeropoint=25,skyrad=[40,50],badpix=[-12000,60000],exact=True)
# 	gauss,psf,psfmag = getpsf.getpsf(data,xpos,ypos,mag,sky,1,1,np.arange(len(xpos)),11,3,'output_psf.fits')
# 	psf_file, hpsf = rdpsf.rdpsf('output_psf.fits')
# 	os.system('mv output_psf.fits %s_psf.fits' %file[:-5])
# 	gauss = [hpsf['GAUSS1'],hpsf['GAUSS2'],hpsf['GAUSS3'],hpsf['GAUSS4'],hpsf['GAUSS5']]
# 	pk = pkfit.pkfit_class(data,gauss,psf_file,1,1)
# 	dummy=0
# 	starflux=[None]*len(xpos)
# 	stardflux=[None]*len(xpos)
# 	for x,y,s in zip(xpos,ypos,sky):
# 		errmag,chi,sharp,niter,scale=pk.pkfit(1,x,y,s,5)
# 		starflux[dummy]=float(scale*10**(0.4*(25.-hpsf['PSFMAG'])))
# 		stardflux[dummy]=float(errmag*10**(0.4*(25.-hpsf['PSFMAG'])))
# 		print('PSF fit to coords %.2f,%.2f gives flux %s +/- %s'%(x,y,starflux[dummy],stardflux[dummy]))
# 		dummy=dummy+1
		
	
	RA=header['RA']
	Dec=header['DEC']
	c=SkyCoord(RA,Dec,frame='icrs',unit=(u.deg,u.deg))
	#Use try and except to look for alternatives to DATE-OBS, like JD,
	try:
		date=header['DATE-OBS']
		t=Time(date,format='isot',scale='utc')
	except:
		try:
			date=header['JD']
			t=Time(date,format='jd',scale='utc')
		except:
			print "Couldn't find header keyword for date.  Quitting...\n"
			quit()
	diff=pyasl.helio_jd(t.jd-2.4e6,c.ra.deg,c.dec.deg,TIME_DIFF=True)
	tdelta=TimeDelta(diff, format='sec')
	hjd=t.jd+tdelta.jd
	
	magtar= -2.5*log10(starflux[0])
	magc1 = -2.5*log10(starflux[1])
	magc2 = -2.5*log10(starflux[2])

	tarc1=magtar-magc1
	tarc2=magtar-magc2

	tarc1err = 2.5*sqrt((stardflux[0]/starflux[0])**2 + (stardflux[1]/starflux[1])**2)
	tarc2err = 2.5*sqrt((stardflux[0]/starflux[0])**2 + (stardflux[2]/starflux[2])**2)

	#now calculate a residual to confirm comparison stars aren't variable
	magcsub = magc1 - magc2
	errmagc = 2.5*sqrt((stardflux[1]/starflux[1])**2 + (stardflux[2]/starflux[2])**2)
	
	print "Finished processing file: %s" %file
	print >>differential, hjd, tarc1, tarc1err, magcsub, errmagc
	
