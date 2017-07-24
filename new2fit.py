import os
import glob
import montage_wrapper as montage

for file in glob.iglob('*.new'):
	name=file[:-4]
	os.system('mv %s.new %snew.fit'%(name,name))
	os.system('rm %s-*'%name)
	os.system('rm %s.solved'%name)
	os.system('rm %s.rdls'%name)
	os.system('rm %s.axy'%name)
	os.system('rm %s.corr'%name)	
	os.system('rm %s.match'%name)	
	os.system('rm %s.wcs'%name)
	
os.system('mkdir astrometry')
os.system('mv *new.fit astrometry/.')
montage.mosaic('astrometry', 'mosaic', background_match=True)