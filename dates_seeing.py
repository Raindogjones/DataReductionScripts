#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import glob
from os.path import exists
from astropy.time import Time
from time import gmtime, strftime
from astropy import units as u
from astropy.coordinates import SkyCoord
from PyAstronomy import pyasl
from astropy.time import TimeDelta
import numpy as np
import photutils
import statistics


dates = open("dates.dat",'w')
for file in glob.iglob('*.fit'):
	hdulist=fits.open(file)
	data=hdulist[0].data
	mean, median, std =sigma_clipped_stats(data,sigma=5.0)
	threshold = mean + 10.5*std
	sources=photutils.irafstarfind(data,threshold=threshold, fwhm=6.0)
	seeing=statistics.median(sources["fwhm"])
	RA=hdulist[0].header['RA']
	Dec=hdulist[0].header['DEC']
	date=hdulist[0].header['DATE-OBS']
	c=SkyCoord(RA,Dec,frame='icrs',unit=(u.deg,u.deg))
	t=Time(date,format='isot',scale='utc')
	diff=pyasl.helio_jd(t.jd-2.4e6,c.ra.deg,c.dec.deg,TIME_DIFF=True)
	tdelta=TimeDelta(diff, format='sec')
	hjd=t.jd+tdelta.jd
	print >> dates, "%s %0.6f %0.2f" %(file,hjd,seeing)