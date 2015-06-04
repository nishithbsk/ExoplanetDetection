from coaddition import *
import os

MAIN_PATH = '/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/eclipse_data/process_output/segments/'
all_segments = os.listdir(MAIN_PATH)
workdir = '/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/eclipse_data/reduction_output/'
ref_cat_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/eclipse_data/reduction_output/ref_catalog.cat'
catalogdir = workdir + 'ref_catalog.cat'

for i, segment in enumerate(all_segments):
	IMAGES_PATH = MAIN_PATH + segment
	setupWorkdir(IMAGES_PATH, workdir)
	detectStars(workdir)
	measureSeeing(workdir)
	createRefCatalog(workdir, outputcatpath=ref_cat_path)
	visualizeAlignmentStars(workdir)
	matchAlignmentStars(workdir, refcat=ref_cat_path)
	createTransforms(workdir)
	transformImages(workdir)
	transformMasks(workdir)
	filename = 'coadd_%d.fits' % i
	filename = workdir + filename
	coaddImages(coaddfile=filename, workdir=workdir)
