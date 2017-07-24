import matplotlib.pyplot as plt
import numpy as np

from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits
import os
from os.path import exists
import sys

import ccdproc

from msumastro import ImageFileCollection, TableTree

nddata.conf.warn_unsupported_correlated = False
imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())

#RED CCD
#gain = 0.8 * u.electron / u.adu
#BLUE CCD
gain = 0.45 * u.electron / u.adu
filt='R_SPECIAL'
object='M3-2'

def oscan_and_trim(image_list):
    """
    Remove overscan and trim a list of images. The original list is replaced by a list of images
    with the changes applied.
    """
    for idx, img in enumerate(image_list):
        oscan = ccdproc.subtract_overscan(img, img[:, 3075:3079], add_keyword={'oscan_sub': True, 'calstat': 'O'}, model=models.Polynomial1D(1))
        image_list[idx] = ccdproc.trim_image(oscan[:, :3073], add_keyword={'trimmed': True, 'calstat': 'OT'})

def trim(image_list):
    """
    Just trim input list of images (bias must be removed by subtraction rather than using overscan). The original list is replaced by a list of images
    with the changes applied.
    """
    for idx, img in enumerate(image_list):
        image_list[idx] = ccdproc.trim_image(img, fits_section='[190:1860,10:960]', add_keyword={'trimmed': True, 'calstat': 'OT'})


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
    
data_dir = '.'

images = ImageFileCollection(data_dir, keywords='*')

if exists ('bias.fit'):
	master_bias=ccdproc.CCDData.read('bias.fit',unit=u.adu)
else:
	bias_list = []
	for hdu, fname in images.hdus(OBJECT='BIAS', return_fname=True):
		meta = hdu.header
		meta['filename'] = fname
		bias_list.append(ccdproc.CCDData(hdu.data,meta=meta, unit=u.adu))
    
	if not bias_list:
		print "Couldn't find any biases. Quitting!"
		sys.exit()

	trim(bias_list)
	biases = ccdproc.Combiner(bias_list)
	master_bias = biases.median_combine(median_func=np.ma.median)
	bias_hdu_list=master_bias.to_hdu()
	bias_hdu_list[0].header["OBJECT"]="MasterCal"
	bias_hdu_list.writeto('bias.fit')

if exists ('%s_flat.fit'%filt):
	master_flat_electron=ccdproc.CCDData.read('%s_flat.fit'%filt,unit=u.electron)
else:
	flats = []
	for flat, fname in images.hdus(OBJECT='FLAT,SKY', return_fname=True):
		meta = flat.header
		meta['filename'] = fname
		if meta["HIERARCH ESO INS FILT1 NAME"]==filt:
			flats.append(ccdproc.CCDData(flat.data, meta=meta, unit=u.adu)) 
    		
	if not flats:
		print "Couldn't find any flats. Quitting!"
		sys.exit()
	

	trim(flats)
	for flat in flats:
		flat_exposure = flat.header['exptime']
		flat = ccdproc.subtract_bias(flat, master_bias,add_keyword={'calib': 'subtracted bias'})
                              

	flat_combiner = ccdproc.Combiner(flats)
	flat_combiner.sigma_clipping(func=med_over_images)
	scaling_func = lambda arr: 1/np.ma.average(arr)
	flat_combiner.scaling = scaling_func
	master_flat = flat_combiner.median_combine(median_func=np.ma.median)
	master_flat.header = flats[0].meta #combiner does not combine metadata!
	master_flat_electron = ccdproc.gain_correct(master_flat, gain=gain)

	flat_hdu_list=master_flat.to_hdu()
	flat_hdu_list[0].header["OBJECT"]="MasterCal"
	flat_hdu_list.writeto('%s_flat.fit'%filt,output_verify='fix')


star_list = []
for star, fname in images.hdus(OBJECT=object, return_fname=True):
    meta = star.header
    meta['filename'] = fname
    if meta["HIERARCH ESO INS FILT1 NAME"]==filt and 'Red' not in fname:
	    star_list.append(ccdproc.CCDData(star.data, meta=meta, unit=u.adu))

if not star_list:
	print "Couldn't find any unreduced object frames. Quitting!"
	sys.exit()

trim(star_list)


star_calibrated = []
for star in star_list:
#	print len(star.data)
#	print star.header['filename']
	fname_base = os.path.basename(star.header['filename'])
	fname_base=fname_base.replace(":","")
	fname_base=fname_base.replace("-","")
	fname_base=fname_base.replace(".","",2)
	fname_base='Red'+fname_base
	if not exists (fname_base):
		star_exp = star.meta['exptime']
		star_bias = ccdproc.subtract_bias(star, master_bias)    
		star_gain = ccdproc.gain_correct(star_bias, gain=gain)
#		star_clean = ccdproc.cosmicray_lacosmic(star_gain, thresh=6, mbox=11, rbox=11, gbox=5)
		star_gain.uncertainty=ccdproc.create_deviation(star_gain,readnoise=3*u.electron).data
#		star_clean = ccdproc.cosmicray_median(star_gain,rbox=5,thresh=0.4)
		new_master_flat=ccdproc.CCDData(master_flat_electron.data, unit=u.electron)
#		star_flat = star_gain.divide(new_master_flat)
		star_flat = ccdproc.flat_correct(star_gain, master_flat_electron)
#		star_flat = ccdproc.flat_correct(star_clean, new_master_flat)
		star_calibrated.append(star_flat)
		image_hdu_list = star_flat.to_hdu()
		image_hdu_list[0].writeto(fname_base,output_verify='fix')
		print "Finished processing %s" %fname_base
	else:
		print "%s already exists. Skipping" %fname_base
    
# for star in star_calibrated:
# 	image_hdu_list = star.to_hdu()
# 	fname_base = os.path.basename(star.header['filename'])
# 	fname_base=fname_base.replace(":","")
# 	fname_base=fname_base.replace("-","")
# 	fname_base=fname_base.replace(".","",2)
# 	image_hdu_list.writeto('Red' + fname_base)
	




