from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import alipy
import glob
import Tkinter as tk
#import ttk
import os
from os.path import exists
from astropy.time import Time
#import astropy.time as time
from astropy import units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
import astropy.coordinates as astrocoords
from astropy.coordinates import EarthLocation
#from PyAstronomy import pyasl
#from time import gmtime, strftime
from astropy.time import TimeDelta
from numpy import log10
from math import sqrt

# from PythonPhot import getpsf
# from PythonPhot import aper
# from PythonPhot import pkfit
# from PythonPhot import rdpsf

#import photutils
# from astropy.stats import sigma_clipped_stats
# from photutils.psf import create_prf
# from photutils.psf import subtract_psf
# from photutils.psf import psf_photometry
# from photutils.psf import GaussianPSF
#from photutils import aperture_photometry, CircularAperture, CircularAnnulus

#import statistics
#from astropy.table.operations import hstack

import sep

import jplephem
import de423



#Things which can be changed go here:
#The aperture radius in arcsec for extraction
arcsec_rad=1.5
selectnewcomp='n' #Change to anything but 'y' in order to use previous comparison selections




#Definitions for accurate calculation of HJD
NTT_LAT = -29.2588
NTT_LON = -70.7339
NTT_ALT = 2375
NTT=EarthLocation.from_geodetic(lat=NTT_LAT, lon=NTT_LON, height=NTT_ALT)

#jd_corr definition from Leo Huckvale's scripts
def jd_corr(mjd, ra, dec, jd_type='bjd'):
     """Return BJD or HJD for input MJD(UTC)."""

     # Initialise ephemeris from jplephem
     eph = jplephem.Ephemeris(de423)

     # Source unit-vector
     ## Assume coordinates in ICRS
     ## Set distance to unit (kilometers)
     src_vec = astrocoords.ICRS(ra=ra, dec=dec,distance=astrocoords.Distance(1, u.km))

     # Convert epochs to astropy.time.Time
     ## Assume MJD(UTC)
     t = Time(mjd, scale='utc', format='mjd', location=NTT)

     # Get Earth-Moon barycenter position
     ## NB: jplephem uses Barycentric Dynamical Time, e.g. JD(TDB)
     ## and gives positions relative to solar system barycenter
     barycenter_earthmoon = eph.position('earthmoon', t.tdb.jd)

     # Get Moon position vectors
     moonvector = eph.position('moon', t.tdb.jd)

     # Compute Earth position vectors
     pos_earth = (barycenter_earthmoon - moonvector * eph.earth_share)*u.km

     if jd_type == 'bjd':
         # Compute BJD correction
         ## Assume source vectors parallel at Earth and Solar System 
		 ## Barycenter
         ## i.e. source is at infinity
         corr = np.dot(pos_earth.T, src_vec.cartesian.xyz)/const.c
     elif jd_type == 'hjd':
         # Compute HJD correction via Sun ephemeris
         pos_sun = eph.position('sun', t.tdb.jd)*u.km
         sun_earth_vec = pos_earth - pos_sun
         corr = np.dot(sun_earth_vec.T, src_vec.cartesian.xyz)/const.c

     # TDB is the appropriate time scale for these ephemerides
     dt = TimeDelta(corr, scale='tdb', format='jd')

     # Compute and return HJD/BJD as astropy.time.Time
     new_jd = t + dt
     return new_jd


#First select central and comparison stars to 

#Defining the pop-up message style
LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Helvetica", 10)
SMALL_FONT = ("Helvetica", 8)
def popupmsg(title,msg):
    popup = tk.Tk()
    popup.wm_title(title)
    label = tk.Label(popup, text=msg, font=LARGE_FONT)
    label.pack(side="top", fill="x", pady=10, padx=10)
    B1 = tk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()



if selectnewcomp=='y':
  

	popupmsg("GAIA", "Use Pick Object function to select target and at \nleast two comparisons, saving their coordinates \neach time to GaiaPick.log \nClick okay to open ref.fit in Gaia and begin selection")
	

	# Open ref.fit to select stars, delete previous star selection every time.
	try:
		os.remove('GaiaPick.Log')
		os.remove('StarPos.dat')
	except:
		pass
		
	os.system('/star/bin/gaia/gaia.sh ref.fit')

xpos,ypos=np.array([0]),np.array([0])
apcoords=list()

gaialog="GaiaPick.Log"
log=open(gaialog)
for line in log:
	if line.startswith('ref.fit'):
		elements=line.split()
		xpos=np.append(xpos,float(elements[1])/2)
		ypos=np.append(ypos,float(elements[2])/2)
		apcoords.append((float(elements[1]),float(elements[2])))

log.close()
xpos=np.delete(xpos,0,0)
ypos=np.delete(ypos,0,0)
xposp=np.zeros(len(xpos))
yposp=np.zeros(len(xpos))

coords=np.stack((xpos,ypos),axis=-1)
coordsp=np.stack((xposp,yposp),axis=-1)


#Now into the image processing
#Find alignment for new images, but first check to make sure that we aren't aligning anything that is already aligned.

images=sorted(glob.glob("Red*.fits"))
try:
	imaligned,x1,y1,x2,y2,x3,y3=np.loadtxt("StarPos.dat",unpack=True)
	imaligned=imaligned.tolist()
except:
	imaligned=[]
	os.system("touch StarPos.dat")

images_to_align = [item for item in images if item not in imaligned]
ref_image="ref.fit"

if not exists (ref_image):
	print "No reference image found for alignment and/or star selection"
	quit()

starpos = open("StarPos.dat",'a')


if images_to_align:
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications: 
		if id.ok == True:
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
			xposp,yposp=id.trans.inverse().apply((xpos,ypos))
			print >>starpos, id.ukn.name, xposp[0], yposp[0], xposp[1], yposp[1], xposp[2], yposp[2]
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
                
starpos.close()
                
#	outputshape = alipy.align.shape(ref_image)

#	for id in identifications:
#		if id.ok == True:
#			alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)

files=np.genfromtxt('StarPos.dat',dtype=None)


#stats = open("stats.dat",'w')
#os.chdir('alipy_out')
#os.system('rm *psf.fits')
differential = open("differential.dat",'w+')
print >>differential, '#HJD \t Tar-c \t Tar-c_err \t comp \t comp_err'

#Parameters for SEP
pixscale=0.24
ganancia=0.91
aper=arcsec_rad/pixscale

aparea=np.pi*aper*aper
ann_in=2.5*aper
ann_out=3.0*aper
ann_area=np.pi*((ann_out*ann_out)-(ann_in*ann_in))



for file in files:
	hdu=fits.open(file[0]+'.fits')
	header=hdu[0].header
	data=hdu[0].data
	x=np.array([file[1],file[3],file[5]])
	y=np.array([file[2],file[4],file[6]])
#Sep commands
	data = data.byteswap().newbyteorder()
	bkg = sep.Background(data)
#	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),xpos+0.5,ypos+0.5,aper,err=bkg.globalrms,gain=ganancia)
#Correction should be -1, I think...
	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),x-1,y-1,aper,err=bkg.globalrms,gain=ganancia)
	nebannflux, nebanndflux, flag=sep.sum_circann(data-bkg.back(),x[0]-1,y[0]-1,ann_in,ann_out,err=bkg.globalrms,gain=ganancia)
# 	starflux,stardflux,flag=sep.sum_circle(data,xpos-1,ypos-1,aper,err=bkg.globalrms,gain=ganancia)
# 	nebannflux, nebanndflux, flag=sep.sum_circann(data,xpos[0]-1,ypos[0]-1,ann_in,ann_out,err=bkg.globalrms,gain=ganancia)
	nebflux=nebannflux*aparea/ann_area
	nebdflux=nebanndflux*aparea/ann_area
	starflux[0]=starflux[0]-nebflux
	stardflux[0]=np.sqrt(stardflux[0]**2 + nebdflux**2)
#	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),xpos-1,ypos-1,aper,err=bkg.globalrms,gain=ganancia,bkgann=(aper+1, aper+2))
#	print file, starflux
	
	
	
	RA=header['RA']
	Dec=header['DEC']
	c=SkyCoord(RA,Dec,frame='icrs',unit=(u.deg,u.deg))
	#Use try and except to look for alternatives to DATE-OBS, like JD,
	try:
		date=header['JD']
		t=Time(date,format='mjd',scale='utc')
#		date=header['DATE-OBS']
#		t=Time(date,format='isot',scale='utc')
	except:
		try:
			date=header['DATE-OBS']
			t=Time(date,format='isot',scale='utc')
#			date=header['JD']
#			t=Time(date,format='jd',scale='utc')
		except:
			print "Couldn't find header keyword for date.  Quitting...\n"
			quit()
	exptime=header['EXPTIME']
	dt=TimeDelta(exptime/2,format='sec')
	tmid=t+dt

	hjd=jd_corr(tmid.mjd, c.ra, c.dec, jd_type='hjd')	
# 	diff=pyasl.helio_jd(tmid.jd-2.4e6,c.ra.deg,c.dec.deg,TIME_DIFF=True)
# 	tdelta=TimeDelta(diff, format='sec')
# 	hjd=tmid.jd+tdelta.jd
	
	magtar= -2.5*log10(starflux[0])
	magc1 = -2.5*log10(starflux[1])
	magc2 = -2.5*log10(starflux[2])

	tarc1=magtar-magc1
	tarc2=magtar-magc2

	tarc1err = 2.5*sqrt((stardflux[0]/starflux[0])**2 + (stardflux[1]/starflux[1])**2)
	tarc2err = 2.5*sqrt((stardflux[0]/starflux[0])**2 + (stardflux[2]/starflux[2])**2)

	print file, hjd, starflux[0], magtar, magc1, tarc1err

	#now calculate a residual to confirm comparison stars aren't variable
	magcsub = magc1 - magc2
	errmagc = 2.5*sqrt((stardflux[1]/starflux[1])**2 + (stardflux[2]/starflux[2])**2)
	
#	print file,magcsub,errmagc
	
	print "Finished processing file: %s" %file
	if tarc1err < 0.2 and not tarc1=="nan" :
		print >>differential, "%.8f %.5f %.5f %.5f %.5f" %(hjd.jd[0], tarc1, tarc1err, magcsub, errmagc)
#	print >>differential, "%.8f %.5f %.5f %.5f %.5f" %(hjd.jd[0], tarc1, tarc1err, magcsub, errmagc)
#	print >>differential, "%.8f %.5f %.5f %.5f %.5f" %(hjd, tarc1, tarc1err, magcsub, errmagc)
	
differential.close()

#Now show a plot! And output the usual differential and comparison plots (must add plot for seeing)
with open ("differential.dat") as dif:
	diff=dif.read()


plotjd,plotmag,ploterr,compmag,comperr=np.array([0]),np.array([0]),np.array([0]),np.array([0]),np.array([0])
diff=diff.split('\n')


for line in diff:
	if not line.startswith('#') and line.startswith('2'):
		line=line.rstrip()
		entries=line.split(' ')
		plotjd=np.append(plotjd,float(entries[0]))
		plotmag=np.append(plotmag,float(entries[1]))
		ploterr=np.append(ploterr,float(entries[2]))		
		compmag=np.append(compmag,float(entries[3]))
		comperr=np.append(comperr,float(entries[4]))		

plotjd=np.delete(plotjd,0,0)
plotmag=np.delete(plotmag,0,0)
ploterr=np.delete(ploterr,0,0)
compmag=np.delete(compmag,0,0)
comperr=np.delete(comperr,0,0)

plt.figure()
plt.errorbar(plotjd,compmag,yerr=comperr,linestyle="None",capsize=2)
plt.gca().invert_yaxis()
plt.ylabel('Magnitude')
plt.xlabel('HJD')
plt.title('Comparison')
plt.savefig('comparison.eps')
plt.clf()

#plt.figure()
plt.errorbar(plotjd,plotmag,yerr=ploterr,linestyle="None",capsize=2)
plt.gca().invert_yaxis()
plt.ylabel('Magnitude')
plt.xlabel('HJD')
plt.title('Differential')
plt.savefig('differential.eps')
plt.show()

