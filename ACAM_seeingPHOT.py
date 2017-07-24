from astropy.io import fits
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import alipy
import glob
import Tkinter as tk
import ttk
import os
from os.path import exists
from astropy.time import Time
#import astropy.time as time
from astropy import units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
import astropy.coordinates as astrocoords
from astropy.coordinates import EarthLocation
from PyAstronomy import pyasl
#from time import gmtime, strftime
from astropy.time import TimeDelta
from math import log10
from math import sqrt

from imexam.imexamine import Imexamine
from imexam import math_helper

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
arcsec_rad=3
seeing_rad=1.5
selectnewcomp='n' #Change to anything but 'y' in order to use previous comparison selections





#Definitions for accurate calculation of HJD
NTT_LAT = -29.2588
NTT_LON = -70.7339
NTT_ALT = 2375
NTT=EarthLocation.from_geodetic(lat=NTT_LAT, lon=NTT_LON, height=NTT_ALT)

INT_LAT = 28.762
INT_LON = -17.878
INT_ALT = 2336
INT=EarthLocation.from_geodetic(lat=INT_LAT, lon=INT_LON, height=INT_ALT)

WHT_LAT = 28.761
WHT_LON = -17.882
WHT_ALT = 2332
WHT=EarthLocation.from_geodetic(lat=WHT_LAT, lon=WHT_LON, height=WHT_ALT)

#jd_corr definition from Leo Huckvale's scripts
def jd_corr(mjd, ra, dec, jd_type='bjd'):
     """Return BJD or HJD for input MJD(UTC)."""

     # Initialise ephemeris from jplephem
     eph = jplephem.Ephemeris(de423)

     # Source unit-vector
     ## Assume coordinates in ICRS
     ## Set distance to unit (kilometers)
#     src_vec = astrocoords.ICRS(ra=ra, dec=dec,unit=(u.degree, u.degree),distance=astrocoords.Distance(1, u.km))
     src_vec = astrocoords.ICRS(ra=ra, dec=dec,distance=astrocoords.Distance(1, u.km))

     # Convert epochs to astropy.time.Time
     ## Assume MJD(UTC)
#     t = Time(mjd, scale='utc', format='mjd', location=NTT)
#     t = Time(mjd, scale='utc', format='mjd', location=WHT)
     t = Time(mjd, scale='utc', format='mjd', location=WHT)

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




#Now into the image processing
#Align new images, but first check to make sure that we aren't aligning anything that is already aligned (and present in the alipy_out folder

images=sorted(glob.glob("RED*.fit"))
images2=sorted(glob.glob("alipy_out/RED*.fits"))
images2=[w.replace('alipy_out/','') for w in images2]
images2=[w.replace('_gregister.fits','.fit') for w in images2]
images_to_align = [item for item in images if item not in images2]
ref_image="ref.fits"

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


#Defining the pop-up message style
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
    
if selectnewcomp=='y':
  
	popupmsg("GAIA", "Use Pick Object function to select target and at \nleast two comparisons, saving their coordinates \neach time to GaiaPick.log \nClick okay to open ref.fits in Gaia and begin selection")

	# Open ref.fit to select stars, delete previous star selection every time.
	try:
		os.remove('GaiaPick.Log')
	except:
		pass
		
	os.system('/star/bin/gaia/gaia.sh ref.fits')

xpos,ypos=np.array([0]),np.array([0])
apcoords=list()

gaialog="GaiaPick.Log"
log=open(gaialog)
for line in log:
	if line.startswith('ref.fits'):
		elements=line.split()
		xpos=np.append(xpos,float(elements[1]))
		ypos=np.append(ypos,float(elements[2]))
		apcoords.append((float(elements[1]),float(elements[2])))

log.close()
xpos=np.delete(xpos,0,0)
ypos=np.delete(ypos,0,0)

coords=np.stack((xpos,ypos),axis=-1)

#stats = open("stats.dat",'w')
os.chdir('alipy_out')
#os.system('rm *psf.fits')
differential = open("differential.dat",'w+')
print >>differential, '#HJD \t Tar-c \t Tar-c_err \t comp \t comp_err \t Seeing \t Airmass'
#print >>differential, '#HJD \t Tar-c \t Tar-c_err \t comp \t comp_err'


pixscale=0.25
ganancia=0.92
aper=arcsec_rad/pixscale


for file in glob.iglob('*.fits'):
	hdu=fits.open(file)
	header=hdu[0].header
	data=hdu[0].data

#commands to get seeing
	test=Imexamine()
	seeing=np.zeros(3)
	xcen=np.zeros(3)
	ycen=np.zeros(3)
	for x in range(0,3):
		monkey=test.gauss_center(int(xpos[x]),int(ypos[x]),data,delta=15)
		seeing[x]=math_helper.gfwhm(monkey[3]) #in pixels
		xcen[x]=monkey[1] #x-centroid in python units
		ycen[x]=monkey[2] #y-centroid in python units
#	test.showplt()

	aper=np.median(abs(seeing))*seeing_rad
#Sep commands
	data = data.byteswap().newbyteorder()
	bkg = sep.Background(data)
	
#Check to make sure centroiding was successful
	if xcen[0]==0 or abs(xcen[0]-xpos[0])>3:
		print "Couldn't centroid %s" %file
		continue
	elif xcen[1]==0 or abs(xcen[1]-xpos[1])>3:
		print "Couldn't centroid %s" %file
		continue
	elif xcen[2]==0 or abs(xcen[2]-xpos[2])>3:
		print "Couldn't centroid %s" %file
		continue

	if ycen[0]==0 or abs(ycen[0]-ypos[0])>3:
		print "Couldn't centroid %s" %file
		continue
	elif ycen[1]==0 or abs(ycen[1]-ypos[1])>3:
		print "Couldn't centroid %s" %file
		continue
	elif ycen[2]==0 or abs(ycen[2]-ypos[2])>3:
		print "Couldn't centroid %s" %file
		continue

	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),xcen-1,ycen-1,aper,err=bkg.globalrms,gain=ganancia)
		
		
#	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),xpos+0.5,ypos+0.5,aper,err=bkg.globalrms,gain=ganancia)
#Correction should be -1, I think...

#	USE THIS ONE FOR FIXED APERTURE!
#	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),xpos-1,ypos-1,aper,err=bkg.globalrms,gain=ganancia)

#	starflux,stardflux,flag=sep.sum_circle(data-bkg.back(),xpos-1,ypos-1,aper,err=bkg.globalrms,gain=ganancia,bkgann=(aper+1, aper+2))
	
	
	RA=header['RA']
	Dec=header['DEC']
	c=SkyCoord(RA,Dec,frame='icrs',unit=(u.deg,u.deg))
	#Use try and except to look for alternatives to DATE-OBS, like JD,
	try:
		date=header['JD']
		t=Time(date,format='jd',scale='utc')
	except:
		try:
			date=header['DATE-OBS']
			t=Time(date,format='isot',scale='utc')
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

	#now calculate a residual to confirm comparison stars aren't variable
	magcsub = magc1 - magc2
	errmagc = 2.5*sqrt((stardflux[1]/starflux[1])**2 + (stardflux[2]/starflux[2])**2)
	
	print "Finished processing file: %s" %file
	if not tarc1=="nan" :
		arcseeing=np.median(abs(seeing))*pixscale
		airmass=header['AIRMASS']
		print >>differential, "%.8f %.5f %.5f %.5f %.5f %.2f %.5f" %(hjd.jd[0], tarc1, tarc1err, magcsub, errmagc,arcseeing,airmass)
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
plt.errorbar(plotjd,compmag,yerr=comperr,linestyle="None")
plt.gca().invert_yaxis()
plt.ylabel('Magnitude')
plt.xlabel('HJD')
plt.title('Comparison')
plt.savefig('comparison.eps')
plt.clf()

#plt.figure()
plt.errorbar(plotjd,plotmag,yerr=ploterr,linestyle="None")
plt.gca().invert_yaxis()
plt.ylabel('Magnitude')
plt.xlabel('HJD')
plt.title('Differential')
plt.savefig('differential.eps')
plt.show()
