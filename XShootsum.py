from astropy.io import fits
import glob
import operator
import numpy as np
from scipy.interpolate import interp1d
import astropy.constants as const
from specutils import Spectrum1D
import astropy.units as u
from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel, Box1DKernel

files = glob.glob("ADP*.fits")
#Assume all files have same wavelength range and set up wavelength and data arrays from that
f=fits.open(files[0])
data=f[1].data
wave=data['WAVE'][0]
flux=data['FLUX'][0]
#Initiate total and running sum to zero
sum=flux*0
err=flux*0
totalerrsq=err
total=sum

nfiles=0
for file in glob.iglob("ADP*.fits"):
	f=fits.open(file)
	data=f[1].data
	header=f[0].header
	wave=data['WAVE'][0]
	flux=data['FLUX'][0]*1e16
	err=data['ERR'][0]*1e16
	total=total+flux
	totalerrsq=totalerrsq+err**2
	nfiles=nfiles+1

average=total/float(nfiles)
averr=np.sqrt(totalerrsq)

outfile=open("average.dat",'w+')
xcor=np.array([wave,average,averr])
xcor=xcor.T
np.savetxt(outfile,xcor,fmt=['%.5f','%.5f','%.5f'])
outfile.close()

