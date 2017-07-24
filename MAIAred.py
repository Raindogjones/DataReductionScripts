#import matplotlib.pyplot as plt
import numpy as np

from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits
import os

import ccdproc

#from msumastro import ImageFileCollection, TableTree

nddata.conf.warn_unsupported_correlated = False
imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())

#Set up only for SLOW readout
gain = 1.25 * u.electron / u.adu

def trim(img,ystart,yend):
	trimmed = ccdproc.trim_image(img, fits_section='[:,%i:%i]'%(ystart,yend))
	return trimmed

def avg_over_images(masked_arr, axis=0):
    """
    Calculate average pixel value along specified axis
    """
    return np.ma.mean(masked_arr, axis=axis)

def med_over_images(masked_arr, axis=0):
    """
    Calculate median pixel value along specified axis
    
    """
    
    dat = masked_arr.data.copy()
    dat[masked_arr.mask] = np.NaN
    return np.ma.median(dat, axis=axis)

#Should be executed in directory with all raw data  
data_dir = '.'

images = ccdproc.ImageFileCollection(data_dir, keywords='*')

# bias_list_1 = []
# bias_list_2 = []
# bias_list_3 = []
# 
# for hdu, fname in images.hdus(IMAGETYP='ZERO', return_fname=True):
# 	meta = hdu.header
# 	meta['filename'] = fname
# 	biasdata1=fits.open(fname)[1].data
# 	bias_list_1.append(ccdproc.CCDData(biasdata1, unit="adu"))
# 	biasdata2=fits.open(fname)[2].data
# 	bias_list_2.append(ccdproc.CCDData(biasdata2, unit="adu"))
# 	biasdata3=fits.open(fname)[3].data
# 	bias_list_3.append(ccdproc.CCDData(biasdata3, unit="adu"))
# 	
# master_bias1 = ccdproc.combine(bias_list_1,method='median')
# master_bias2 = ccdproc.combine(bias_list_2,method='median')
# master_bias3 = ccdproc.combine(bias_list_3,method='median')


#Alternative routine to fudge biases
# bias=np.full((3100,2150),1050,dtype=np.float)
# master_bias1=ccdproc.CCDData(bias,unit="adu")
# master_bias2=ccdproc.CCDData(bias,unit="adu")
# master_bias3=ccdproc.CCDData(bias,unit="adu")
# meta=fits.PrimaryHDU().header
# meta['INSTRUME']="MAIA"
# meta['OBSMODE']="BIAS"
# meta['IMAGETYP']="MasterCal"
# meta['OBJECT']="MasterBias"
# meta['EXPTIME']=0
# meta['EXPTYPE']="BIAS"
# meta['DETMODE']="LGN"
# meta['READMODE']="L"
# meta['DETGAIN']=1.25
# meta['DETBIAS']=1050
# meta['DETNAME']="EDDINGTON-13"
# meta['DETTYPE']="E2V42-C0"
# meta['DETID']="DET99"


# hdu1=fits.ImageHDU(master_bias1)
# hdu2=fits.ImageHDU(master_bias2)
# hdu3=fits.ImageHDU(master_bias3)
# primary=fits.ImageHDU(data=None, header=meta)
# primary.header['SIMPLE']="T"
# primary.header['EXTEND']="T"
# hdulist=fits.HDUList([hdu1,hdu2,hdu3])
# hdulist.verify('fix')
# #hdulist.insert(0,primary)
# hdulist[0].header=meta
# hdulist[0].header["OBJECT"]="MasterBias"
# hdulist[0].header["IMAGETYP"]="MasterCal"
# hdulist.writeto('MasterBias.fits',output_verify='fix',clobber=True)
# 
# if os.path.isfile('MasterBias.fits'):
# 	print "Successfully produced MasterBias frame!\n"

flat_list_1 = []
flat_list_2 = []
flat_list_3 = []
for hdu, fname in images.hdus(IMAGETYP="FLAT", return_fname=True):
	flatmeta = hdu.header
	flatmeta['filename'] = fname
	flatdata1=fits.open(fname)[1].data
	if np.mean(flatdata1) > 10000:
		flat_list_1.append(ccdproc.CCDData(flatdata1, unit="adu"))
	flatdata2=fits.open(fname)[2].data
	if np.mean(flatdata2) > 10000:
		flat_list_2.append(ccdproc.CCDData(flatdata2, unit="adu"))
	flatdata3=fits.open(fname)[3].data
	if np.mean(flatdata3) > 10000:
		flat_list_3.append(ccdproc.CCDData(flatdata3, unit="adu"))


#If using masterbias:
#for flat in flat_list_1:
#	flat=ccdproc.subtract_bias(flat, master_bias1, add_keyword={'calib': 'subtracted_bias'})
#for flat in flat_list_2:
#	flat=ccdproc.subtract_bias(flat, master_bias2, add_keyword={'calib': 'subtracted_bias'})
#for flat in flat_list_3:
#	flat=ccdproc.subtract_bias(flat, master_bias3, add_keyword={'calib': 'subtracted_bias'})

#Using overscan:
for idx, flat in enumerate(flat_list_1):
	flat_list_1[idx]=ccdproc.subtract_overscan(flat, fits_section='[2110:2140,:]')
	flat_list_1[idx]=ccdproc.trim_image(flat, fits_section='[52:2100,:]')
for idx, flat in enumerate(flat_list_2):
	flat_list_2[idx]=ccdproc.subtract_overscan(flat, fits_section='[2110:2140,:]')
	flat_list_2[idx]=ccdproc.trim_image(flat, fits_section='[52:2100,:]')
for idx, flat in enumerate(flat_list_3):
	flat_list_3[idx]=ccdproc.subtract_overscan(flat, fits_section='[2110:2140,:]')
	flat_list_3[idx]=ccdproc.trim_image(flat, fits_section='[52:2100,:]')	

                           
flat1_combiner = ccdproc.Combiner(flat_list_1)
flat1_combiner.sigma_clipping(func=med_over_images)
scaling_func = lambda arr: 1/np.ma.average(arr)
flat1_combiner.scaling = scaling_func
master_flat1 = flat1_combiner.median_combine(median_func=np.ma.median)
master_flat1_electron = ccdproc.gain_correct(master_flat1, gain=gain)

flat2_combiner = ccdproc.Combiner(flat_list_2)
flat2_combiner.sigma_clipping(func=med_over_images)
scaling_func = lambda arr: 1/np.ma.average(arr)
flat2_combiner.scaling = scaling_func
master_flat2 = flat2_combiner.median_combine(median_func=np.ma.median)
master_flat2_electron = ccdproc.gain_correct(master_flat2, gain=gain)

flat3_combiner = ccdproc.Combiner(flat_list_3)
flat3_combiner.sigma_clipping(func=med_over_images)
scaling_func = lambda arr: 1/np.ma.average(arr)
flat3_combiner.scaling = scaling_func
master_flat3 = flat3_combiner.median_combine(median_func=np.ma.median)
master_flat3_electron = ccdproc.gain_correct(master_flat3, gain=gain)
	
hdu1=fits.ImageHDU(master_flat1)
hdu2=fits.ImageHDU(master_flat2)
hdu3=fits.ImageHDU(master_flat3)
hdulist=fits.HDUList([hdu1,hdu2,hdu3])
hdulist.verify('fix')
hdulist[0].header=flatmeta
hdulist[0].header["OBJECT"]="MasterFlat"
hdulist[0].header["IMAGETYP"]="MasterCal"
hdulist.writeto('MasterFlat.fits',output_verify='fix',clobber=True)
if os.path.isfile('MasterFlat.fits'):
	print "Successfully produced MasterFlat frame!\n"


for star, fname in images.hdus(IMAGETYP='OBJECT', return_fname=True):
	meta=star.header
	meta['filename']=fname
	try:
		ys=int(meta['YSTART'][1:-1])
		ye=ys+int(meta['WINY'])-1
		new_master_flat1=trim(ccdproc.CCDData(master_flat1_electron.data,unit=u.electron), ys, ye)
		new_master_flat2=trim(ccdproc.CCDData(master_flat2_electron.data,unit=u.electron), ys, ye)
		new_master_flat3=trim(ccdproc.CCDData(master_flat3_electron.data,unit=u.electron), ys, ye)
	except:
		new_master_flat1=ccdproc.CCDData(master_flat1_electron.data,unit=u.electron)
		new_master_flat2=ccdproc.CCDData(master_flat2_electron.data,unit=u.electron)
		new_master_flat3=ccdproc.CCDData(master_flat3_electron.data,unit=u.electron)
	star1=ccdproc.CCDData(fits.open(fname)[1].data,unit="adu")
	star1_bias=ccdproc.subtract_overscan(star1, fits_section='[2110:2140,:]')
	star1_bias=ccdproc.trim_image(star1_bias, fits_section='[52:2100,:]')
#	star1_bias=ccdproc.subtract_bias(star1, ccdproc.CCDData(trim(master_bias1, ys, ye),unit="adu"))
	star1_gain=ccdproc.gain_correct(star1_bias,gain=gain)
	star1_flat=ccdproc.flat_correct(star1_gain,new_master_flat1)

	head1=fits.open(fname)[1].header
	
	star2=ccdproc.CCDData(fits.open(fname)[2].data,unit="adu")
	star2_bias=ccdproc.subtract_overscan(star2, fits_section='[2110:2140,:]')
	star2_bias=ccdproc.trim_image(star2_bias, fits_section='[52:2100,:]')
#	star2_bias=ccdproc.subtract_bias(star2, ccdproc.CCDData(trim(master_bias2, ys, ye),unit="adu"))
	star2_gain=ccdproc.gain_correct(star2_bias,gain=gain)
	star2_flat=ccdproc.flat_correct(star2_gain,new_master_flat2)
	head2=fits.open(fname)[2].header

	star3=ccdproc.CCDData(fits.open(fname)[3].data,unit="adu")
	star3_bias=ccdproc.subtract_overscan(star3, fits_section='[2110:2140,:]')
	star3_bias=ccdproc.trim_image(star3_bias, fits_section='[52:2100,:]')
#	star3_bias=ccdproc.subtract_bias(star3, ccdproc.CCDData(trim(master_bias3, ys, ye),unit="adu"))
	star3_gain=ccdproc.gain_correct(star3_bias,gain=gain)
	star3_flat=ccdproc.flat_correct(star3_gain,new_master_flat3)
	head3=fits.open(fname)[3].header
	
# 	print "Bias %d" %np.mean(star3_bias)
# 	print "Gain %d" %np.mean(star3_gain)
# 	print "New Flat %d" %np.mean(new_master_flat3)	
# 	print "Flattened %d" %np.mean(star3_flat)
	
	hdu1=fits.ImageHDU(star1_flat)
	hdu2=fits.ImageHDU(star2_flat)
	hdu3=fits.ImageHDU(star3_flat)
	hdulist=fits.HDUList([hdu1,hdu2,hdu3])
	hdulist.verify('fix')
	hdulist[0].header=meta
	hdulist[0].header["IMAGETYP"]="REDOBJECT"
	hdulist[1].header=head1
	hdulist[2].header=head2
	hdulist[3].header=head3
	for x in range(1,4):
		hdulist[x].header['BZERO']=0
	fname_base = os.path.basename(meta['filename'])
	hdulist.writeto('RED' + fname_base,output_verify='fix',clobber=True)
	print "%s reduced and written to file" %meta['filename'][:]


	
# star_list = []
# for star, fname in images.hdus(IMAGETYP='OBJECT', return_fname=True):
# 	meta = star.header
# 	meta['filename'] = fname
# 	if meta["WFFBAND"]==band:
# 		stardata=fits.open(fname)[1].data
# 		star_list.append(ccdproc.CCDData(stardata, meta=meta, unit="adu"))
# 
# 
# 
# star_calibrated = []
# for star in star_list:
# 	star_exp = star.meta['exptime']
# 	star_bias = ccdproc.subtract_bias(star, master_bias)    
# 	star_gain = ccdproc.gain_correct(star_bias, gain=gain)
# #		star_clean = ccdproc.cosmicray_lacosmic(star_gain, thresh=5.5, mbox=11, rbox=11, gbox=5)
# #		star_flat = ccdproc.flat_correct(star_clean, master_flat_electron)
# #		star_flat = ccdproc.flat_correct(star_gain, master_flat_electron)
# #		star_calibrated.append(star_flat)
# 	new_master_flat=ccdproc.CCDData(master_flat_electron.data, unit=u.electron)
# 	star_flat = ccdproc.flat_correct(star_gain, new_master_flat)
# 	image_hdu_list=star_flat.to_hdu()
# 	fname_base = os.path.basename(star.meta['filename'])
# 	image_hdu_list.writeto('RED' + fname_base,output_verify='fix',clobber=True)
# 	print "%s reduced and written to file" %star.meta['filename'][:]






		
		
		
		
		
		