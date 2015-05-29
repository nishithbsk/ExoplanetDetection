from coaddition import *

IMAGES_PATH = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/red_images.list'
workdir = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/red_data_reduction_workdir/'
ref_cat_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/red_data_reduction_workdir/ref_catalog.cat'
catalogdir = workdir + 'ref_catalog.cat'

setupWorkdir(IMAGES_PATH, workdir)
detectStars(workdir)
measureSeeing(workdir)
createRefCatalog(workdir, outputcatpath=ref_cat_path)
visualizeAlignmentStars(workdir)
matchAlignmentStars(workdir, refcat=ref_cat_path)
createTransforms(workdir)
transformImages(workdir)
transformMasks(workdir)
coaddImages(workdir

IMAGES_PATH = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/blue_images.list'
workdir = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/blue_data_reduction_workdir/'
ref_cat_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/red_data_reduction_workdir/ref_catalog.cat'
catalogdir = workdir + 'ref_catalog.cat'

setupWorkdir(IMAGES_PATH, workdir)
detectStars(workdir)
measureSeeing(workdir)
createRefCatalog(workdir, outputcatpath=catalogdir)
visualizeAlignmentStars(workdir)
matchAlignmentStars(workdir, refcat=ref_cat_path)
createTransforms(workdir)
transformImages(workdir)
transformMasks(workdir)
coaddImages(workdir)

IMAGES_PATH = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/green_images.list'
workdir = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/green_data_reduction_workdir/'
ref_cat_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/calibrated_stars/BD/output_BD/red_data_reduction_workdir/ref_catalog.cat'
catalogdir = workdir + 'ref_catalog.cat'

setupWorkdir(IMAGES_PATH, workdir)
detectStars(workdir)
measureSeeing(workdir)
createRefCatalog(workdir, outputcatpath=catalogdir)
visualizeAlignmentStars(workdir)
matchAlignmentStars(workdir, refcat=ref_cat_path)
createTransforms(workdir)
transformImages(workdir)
transformMasks(workdir)
coaddImages(workdir)