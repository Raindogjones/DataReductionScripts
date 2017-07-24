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
gain = 1.6 * u.electron / u.adu

def trim(image_list):
    """
    Just trim input list of images (bias must be removed by subtraction rather than using overscan). The original list is replaced by a list of images
    with the changes applied.
    """
    for idx, img in enumerate(image_list):
#        image_list[idx] = ccdproc.trim_image(img, fits_section='[60:2080,20:4060]', add_keyword={'trimmed': True, 'calstat': 'OT'})
        image_list[idx] = ccdproc.trim_image(img, fits_section='[100:2040,60:4020]', add_keyword={'trimmed': True, 'calstat': 'OT'})


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

bias_list = []
for hdu, fname in images.hdus(IMAGETYP='zero', return_fname=True):
    meta = hdu.header
    meta['filename'] = fname
    if meta['CCDSPEED']=="SLOW":
		biasdata=fits.open(fname)[1].data
		bias_list.append(ccdproc.CCDData(biasdata, unit="adu"))


#trim(bias_list)
master_bias = ccdproc.combine(bias_list,method='median')

bias_hdu_list=master_bias.to_hdu()
bias_hdu_list[0].header["OBJECT"]="MasterCal"
bias_hdu_list[0].header["IMAGETYP"]="MasterBias"
bias_hdu_list.writeto('bias.fits',output_verify='fix',clobber=True)

if os.path.isfile('bias.fits'):
	print "Successfully produced bias frame!\n"

flatbands = list()
for flat, fname in images.hdus(IMAGETYP='sky', return_fname=True):
	meta = flat.header
	if meta["WFFBAND"] not in flatbands:
		flatbands.append(meta["WFFBAND"])
		print "Flats found for filter: %s\n" %meta["WFFBAND"]

for band in flatbands:
	flats=[]
	for flat, fname in images.hdus(IMAGETYP='sky', return_fname=True):
		meta=flat.header
		if meta["WFFBAND"]==band:
			flatdata=fits.open(fname)[1].data
			flats.append(ccdproc.CCDData(flatdata, meta=meta, unit="adu"))
	
#	trim(flats)
	for flat in flats:
		flat_exposure = flat.header['EXPTIME']
    	flat = ccdproc.subtract_bias(flat, master_bias,
                              add_keyword={'calib': 'subtracted bias'})
                              
	flat_combiner = ccdproc.Combiner(flats)
	flat_combiner.sigma_clipping(func=med_over_images)
	scaling_func = lambda arr: 1/np.ma.average(arr)
	flat_combiner.scaling = scaling_func
	master_flat = flat_combiner.median_combine(median_func=np.ma.median)
#	master_flat=ccdproc.combine(flats,method='median')
	master_flat.header = flats[0].meta #combiner does not combine metadata!

	master_flat_electron = ccdproc.gain_correct(master_flat, gain=gain)
	
	flat_hdu_list=master_flat.to_hdu()
	flat_hdu_list[0].header["OBJECT"]="MasterCal"
	flat_hdu_list[0].header["IMAGETYP"]="MasterFlat"
#	flat_hdu_list[0].header["WFFBAND"]=band
	flat_hdu_list.writeto('%s_flat.fits'%band,output_verify='fix',clobber=True)
	if os.path.isfile('%s_flat.fits'%band):
		print "Successfully produced flat for filter %s!\n" %band

	
	star_list = []
	for star, fname in images.hdus(IMAGETYP='object', return_fname=True):
		meta = star.header
		meta['filename'] = fname
		if meta["WFFBAND"]==band:
			stardata=fits.open(fname)[1].data
			star_list.append(ccdproc.CCDData(stardata, meta=meta, unit="adu"))

#	trim(star_list)
	
# 	star_calibrated = []
# 	for star in star_list:
#     	star_exp = star.meta['exptime']
#     	star_bias = ccdproc.subtract_bias(star, master_bias)    
#     	star_gain = ccdproc.gain_correct(star_bias, gain=gain)
#     	star_clean = ccdproc.cosmicray_lacosmic(star_gain, thresh=5, mbox=11, rbox=11, gbox=5)
#     	star_flat = ccdproc.flat_correct(star_clean, master_flat_electron)
#     	star_calibrated.append(star_flat)
# 
#     
# 	for star in star_calibrated:
# 		image_hdu_list = star.to_hdu()
# 		fname_base = os.path.basename(star.header['filename'])
# 		image_hdu_list.writeto('Red' + fname_base)

	star_calibrated = []
	for star in star_list:
		star_exp = star.meta['exptime']
		star_bias = ccdproc.subtract_bias(star, master_bias)    
		star_gain = ccdproc.gain_correct(star_bias, gain=gain)
#		star_clean = ccdproc.cosmicray_lacosmic(star_gain, thresh=5.5, mbox=11, rbox=11, gbox=5)
#		star_flat = ccdproc.flat_correct(star_clean, master_flat_electron)
#		star_flat = ccdproc.flat_correct(star_gain, master_flat_electron)
#		star_calibrated.append(star_flat)
		new_master_flat=ccdproc.CCDData(master_flat_electron.data, unit=u.electron)
		star_flat = ccdproc.flat_correct(star_gain, new_master_flat)
		image_hdu_list=star_flat.to_hdu()
		fname_base = os.path.basename(star.meta['filename'])
		image_hdu_list.writeto('RED' + fname_base,output_verify='fix',clobber=True)
		print "%s reduced and written to file" %star.meta['filename'][:]


# masterflats=[]
# for flat, fname in images.hdus(IMAGETYP='MasterFlat', OBJECT='MasterCal', return_fname=True):
#     meta = flat.header
#     if meta["WFFBAND"] not in flatbands:
#     	masterflats.append(ccdproc.CCDData(flat.data, meta=meta, unit="adu"))
#     	print "Master flat found for filter: %s\n" %meta["WFFBAND"]
#     	
# for flat in masterflats:
# 	star_list = []
# 	for star, fname in images.hdus(IMAGETYP='object', return_fname=True):
#     	meta = star.header
#     	meta['filename'] = fname
#     	if meta["WFFBAND"]==flat.header["WFFBAND"]:
# 	    	star_list.append(ccdproc.CCDData(star.data, meta=meta, unit="adu"))
# 
# 	trim(star_list)
# 	
# 	star_calibrated = []
# 	for star in star_list:
#     	star_exp = star.meta['exptime']
#     	star_bias = ccdproc.subtract_bias(star, master_bias)    
#     	star_gain = ccdproc.gain_correct(star_bias, gain=gain)
#     	star_clean = ccdproc.cosmicray_lacosmic(star_gain, thresh=5, mbox=11, rbox=11, gbox=5)
#     	star_flat = ccdproc.flat_correct(star_clean, flat)
#     	star_calibrated.append(star_flat)
# 
#     
# 	for star in star_calibrated:
# 		image_hdu_list = star.to_hdu()
# 		fname_base = os.path.basename(star.header['filename'])
# 		image_hdu_list.writeto('Red' + fname_base)



		
		
		
		
		
		