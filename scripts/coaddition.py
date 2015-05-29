'''
coaddition.py

This script manages the coaddition of astronomical data. Under the hood,
this script calls the IRAF library to align and resample images. The
status of the coaddition is persistant, so it may be completed over
multiple sessions.

NOTE: If some images are rejected because they cannot be aligned with the
reference, the user should restart the process from the beginning with
those files removed.

IRAF: http://iraf.net

This module is based on ALIPY by Malte Tewes at the University of Geneva,
Switzerland.
http://obswww.unige.ch/~tewes/alipy/
'''

##################################################

import sys, os, shutil, math
import numpy as np
from numpy import *
from pyraf import iraf
from glob import glob
from modules.variousfct import *
from modules.star import *
from datetime import datetime, timedelta
sys.path.append("/afs/ir.stanford.edu/class/physics100/lib/python/modules")
import f2n, pyfits


##################################################

sextractor="/afs/ir/class/physics100/src/ldacpipeline-0.12.20/bin/Linux_64/sex_theli"
alipy_sex_configdir = '/afs/ir.stanford.edu/class/physics100/src/alipy_1.3/sextractor_config'

###################################################
# No Need to call this method in your script
class UnrecognizedArgumentFormatException(Exception): pass

def checkForTextFile(fileargument):
    '''
    checkForTextFile - This is a utility function that allows us to pass
         either a list of filenames or a textfile with one filename per line
         into any of the above functions
         
     @param fileargument - the argument passed to the calling function to
               be inspected
     @returns list of filenames
     '''

    if isinstance(fileargument, str):
        return [ filename.strip() \
                     for filename in open(fileargument).readlines()]
    if isinstance(fileargument, type([])) and isinstance(fileargument[0], str):
        return fileargument
    raise UnrecognizedArgumentFormatException(fileargument)

########################################################################

def setupWorkdir(inputfiles, workdir, forceDeleteWorkdir = False):
    '''
    Does initial organization to keep track of workproducts produced
    during alignment and coaddition. WARNING: WILL REMOVE WORKDIR

    @parameter inputfiles : a textfilename or a list of filepaths specifying
    the images to align and coadd.
    @parameter workdir : the directory where workfiles will be stored,
                  e.g. ~/physics100/workdir/TAgroup/M53/G/coadd

    '''


    rawimgpaths = checkForTextFile(inputfiles)

    print "I have found %i images to align." % len(rawimgpaths)

    # we build a little database (list of dicts) for our images. Will be completed later.
    db = [{'i':i, 'scalingfactor':1.0, 'rawimgpath':rawimgpath, 'rawimgname':os.path.split(rawimgpath)[-1]} for i, rawimgpath in enumerate(rawimgpaths)]

    # we check for the working directory :
    if os.path.isdir(workdir):
        if not forceDeleteWorkdir:
            print "Your workdir already exists. "
            print "(workdir = %s)" % workdir
            toremove = raw_input('Should I delete it?').lower()
            if toremove == 'y' or toremove == 'yes':
                shutil.rmtree(workdir)
                print "Ok, starting from scratch."
                os.mkdir(workdir)
        else:
            print "Your workdir already exists. I will remove it. "
            print "(workdir = %s)" % workdir
            shutil.rmtree(workdir)
            print "Ok, starting from scratch."
            os.mkdir(workdir)
            
    else:
        os.mkdir(workdir)




    # We pickle the database, and proceed with the next script.
    writepickle(db, os.path.join(workdir, "db.pkl"))

    print "Database written. See you !"


####################################################################

def detectStars(workdir):
    '''
    Runs Sextractor to find stars for alignment in each image.

    @parameter workdir : directory where workfiles are stored for
    alignment
    '''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))


    print "I will now start to run sextractor on these %i images." % len(db)

    starttime = datetime.now()

    # for this step we work in the sextractor_config dir :

    if not os.path.isdir(alipy_sex_configdir):
            raise mterror("There is no sextractor_config directory !")
    curdir=os.getcwd()
    os.chdir(alipy_sex_configdir)

    for image in db:

            print "- " * 40
            print image['i']+1, "/", len(db), ":", image['rawimgname']

            # we change the name of the catalog :
            image['basename'] = os.path.splitext(image['rawimgname'])[0] # this is the filename without extension.
            image['catfilename'] = image['basename'] + ".cat"
            image['catfilepath'] = os.path.join(workdir, image['catfilename'])

            sexout = os.system('%s %s -CATALOG_NAME %s' % (sextractor, image['rawimgpath'], image['catfilepath']))

            print image['catfilename']


            # move check image in case it was created
            if os.path.isfile("check.fits"):
                    shutil.move("check.fits", os.path.join(workdir, image['basename'] + "_checkimg.fits"))

    print "- " * 40
    os.chdir(curdir)

    sexdonetime = datetime.now()
    timetaken = nicetimediff(sexdonetime - starttime)
    print "I'me done. It took me %s" % timetaken

    # We pickle the database, and proceed with the next script.
    writepickle(db, os.path.join(workdir, "db.pkl"))

    print "Database written. See you !"


##################################################################

def measureSeeing(workdir, showseeingplots=False):
    '''
    Measures the typical FWHM of stars in each image. This step is needed
    to determine to which image we should align.

    @parameter workdir : directory where workfiles are stored for
                          alignment
    @parameter showseeingplots : If true, will show interactive plots with
                          histograms of star FWHMs, one per
                          image. CURRENTLY BROKEN.
    '''

    
    if showseeingplots:
	print "Importing matplotlib..."
	import matplotlib.pyplot as plt 	

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))


    # Let's show some stats on these catalogs, might be useful...

    for image in db:

            print "- " * 40
            print image['i']+1, "/", len(db), ":", image['rawimgname']

            goodsexstars = readsexcatasstars(image['catfilepath'])
            image['nbrsources'] = len(goodsexstars)
            sortedsexstars = sortstarlistby(goodsexstars, 'fwhm')

            # searching for the mode of the distribution :

            fwhms = np.array([star.fwhm for star in sortedsexstars])
            if len(fwhms) > 10:
                    # In pixels, we expect the fwhm to lie between 1.5 and 7.0 ...
                    nbins = 40
                    (hist, edges) = np.histogram(fwhms, bins=nbins, range=(0, 20))
                    # We find the peak, and build a narrower hist around it
                    maxpos = argmax(hist)
                    print maxpos
                    if maxpos == 0 or maxpos == nbins-1:
                            raise mterror("The FWHM distribution is anormal. Something is wrong with sextractor...")
                    peakpos = 0.5*(edges[maxpos] + edges[maxpos+1])
                    (hist, edges) = np.histogram(fwhms, bins=10, range=(peakpos-1, peakpos+1))
                    maxpos = argmax(hist)
                    seeingpixels = edges[maxpos]+0.1

            elif len(fwhms) > 0: # not many stars ... we use all of them
                    seeingpixels = np.median(fwhms)
            else :
                    seeingpixels = -1.0	

            print "Stellar FWHM [pixels] :", seeingpixels
            image['seeingpixels'] = seeingpixels

            # Optional plot, for tests
            if showseeingplots:
                    plt.figure(image['i']) 
                    plt.hist(fwhms, range=(0, 10), bins=30, facecolor='g', alpha=0.75)
                    plt.xlabel('FWHM [pixels]')
                    plt.title('%s'%image['rawimgname'])
                    plt.axvline(x=seeingpixels, linewidth=2, color='r')
                    plt.figtext(0.5, 0.8, r'$\mathrm{Measured\ FWHM\ [pixels]\ :}\ %5.2f$'%seeingpixels)
                    plt.grid(True)



            # And we measure the ellipticity of the images, by looking at sources with similar width then our seeingpixels

            ells = array([star.ell for star in sortedsexstars if abs(star.fwhm - seeingpixels)/seeingpixels < 0.1])
            print "I found", len(ells), "stars for ellipticity measure."

            if len(ells) > 0:
                    ell = median(ells)	
            else:
                    print "Bummer !"
                    ell = -1.0

            print "Measured ellipticity :", ell
            image['ell'] = ell
    print "- " * 40


    print "\nSo, here are some stats :"
    print " -   full number of good sextractor sources (not necessarily stars),"
    print " -   stellar FHWM in pixels,"
    print " -   median of sextractor's ELLIPTICITY for sources with stellar FWHM.\n"

    for image in db:
            print "| %4i | %6.3f | %6.3f | %s " % (image['nbrsources'], image['seeingpixels'], image['ell'], image['rawimgname'])

    # We "update" the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))

    print "\nYou could use this to choose a reference image for the alignment ..."

    # Just for info, we find one image with the most nbrsources detected by sextractor :
    sorteddb = sorted(db, key=lambda k: k['nbrsources'])
    print "For your information, if you choose autorefselect = True, I would use :"
    print sorteddb[-1]['rawimgname']
    print "as reference image."

    if showseeingplots:
        print "(Type show() to see the plots)"



####################################################################

def createRefCatalog(workdir, outputcatpath = None, 
                                 handpickedcat = None, 
                                 minfwhm = 1.0, maxfwhm =
                                 7.5, maxell = 1.0):
# 0.7
    '''Creates a catalog of the best stars for alignment from the image
    with the best seeing, as identified in measureSeeing
    
    @parameter workdir : directory where workfiles are stored for
                          alignment
    @parameter outputcatpath : path to file where a copy of the reference
                          catalog should be stored. Useful when other filters need to be aligned
                          this this observaiton
    @parameter handpickedcat : name of a txt file with hand selected stars
                          for alignment. Setting this variable overrides
                          automatic star selection.
    @parameter minfwhm : minimum FWHM (in pixels) for a star to be
                          selected as a reference star
    @parameter maxfwhm : maximum FWHM (in pixels) for a star to be
                          selected as a reference star
    @parameter maxell : maximum ellipticity of a star to be selected as a
                          reference star
    '''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    print "Ok, we have %i images." % len(db)

    # Check a few things about ref image :


    print "Automatic reference image selection"
    sorteddb = sorted(db, key=lambda k: k['nbrsources'])
    refimg = sorteddb[-1]
    refimgname = refimg['rawimgname']

    print "Reference image catalog : %s" % refimg['catfilename']


    # Manual or automatic alignment stars selection :

    autofindalistars = handpickedcat is None

    maxnbralistars = 60

    if autofindalistars :
            print "So you want me to select alignment stars. Fine, let's see what we can do."

            goodsources = readsexcatasstars(refimg['catfilepath'])
            print "Parameters minfwhm = %f, maxfwhm = %f, maxell = %f" % (minfwhm, maxfwhm, maxell)
            autostars = [star for star in goodsources if star.fwhm > minfwhm and star.fwhm < maxfwhm and star.ell < maxell]
            autostars = sortstarlistbyflux(autostars)
            autostarssel = autostars[:maxnbralistars]

            print "I've selected the %i brightest stars among %i available good stars." %(len(autostarssel), len(autostars))
            for star in autostarssel:
                    star.write()

            preciserefmanstars = autostarssel

    else :
            manalistarscat = handpickedcat
            print "I will use your hand-picked alignment stars."
            if not os.path.isfile(manalistarscat):
                    raise mterror("Catalog %s does not exist." % manalistarscat)
            refmanstars = readmancatasstars(manalistarscat)

            print "Ok, I will see if I can find these %i stars in the sextractor catalog." % len(refmanstars)

            refautostars = readsexcatasstars(refimg['catfilepath'])
            refautostars = sortstarlistbyflux(refautostars)
            print "The sextractor catalog contains %i usable stars." % len(refautostars)

            print "Identification tolerance : %f pixels." % identificationtolerance
            preciserefmanstars = listidentify(refmanstars, refautostars, identificationtolerance)
            preciserefmanstars = sortstarlistbyflux(preciserefmanstars)

            if len(preciserefmanstars) != len(refmanstars):
                    raise mterror("Could not identify all the alignment stars in the sextractor catalog of ref image.\n(See above to identify problematic stars.)")
            else:
                    print "No problematic stars in your hand-picked selection, good job for a human."


    # We write preciserefmanstars into a pickle, it's the only thing that will be used by the next scripts.
    writepickle(preciserefmanstars, os.path.join(workdir, "alistars.pkl"))
    
    if outputcatpath is not None:
        shutil.copyfile(os.path.join(workdir, "alistars.pkl"), outputcatpath)

    print "So far so good. Reference catalog created."

##########################################################################


def visualizeAlignmentStars(workdir):
    '''Creates an image with the alignment reference stars highlighted,
    for inspection purposes.

    @parameter workdir : directory where workfiles are stored for
                          alignment

    '''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    print "Ok, we have %i images." % len(db)
    print "Reading alignment star catalog..."
    alistars = readpickle(os.path.join(workdir, "alistars.pkl"))

    alistarsdicts = [{"name":star.name, "x":star.x, "y":star.y} for star in alistars]

    # Get the ref image (similar to previous script)
    print "Automatic reference image selection :"
    sorteddb = sorted(db, key=lambda k: k['nbrsources'])
    refimg = sorteddb[-1]
    refimgname = refimg['rawimgname']
    print refimgname


    pnginfostring = "Alignment star selection (%i stars)" % len(alistarsdicts)

    # We are ready to proceed with the png creation

    f2nimage = f2n.fromfits(refimg['rawimgpath'], hdu=0, verbose=True)
    f2nimage.setzscale("auto", "auto")
    f2nimage.rebin(3)


    f2nimage.makepilimage(scale="log", negative=False)
    #f2nimage.drawstarslist(sexgoodsources, r=20)
    f2nimage.drawstarslist(alistarsdicts, r=20)
    f2nimage.writeinfo([pnginfostring])
    f2nimage.writetitle(os.path.basename("Ref : " + refimgname))

    pngfilepath = os.path.join(workdir, "alistars_refimage.png")
    f2nimage.tonet(pngfilepath)

    print "Done."
    print "Now have a look at :"
    print pngfilepath
    print "... and decide if you are happy with this selection."

	

######################################################

def matchAlignmentStars(workdir, refcat = None, tolerance = 5.0,
                        minnbrstars = 10, maxdist = 15):

    '''Matches stars from each image to a reference catalog. Calls IRAF's
    geomap to create a transformation between the original image and the
    reference catalog.

    @parameter workdir : directory where workfiles are stored for
                          alignment
    @parameter refcat : Reference catalog produced by createRefCatalog. If
                          None, defaults to reference catalog in
                          workdir. Set to a /path/to/file if you wish to
                          align these images to an external image, ie a
                          another filter.
    @parameter tolerance : Error allowed (in pixels) when matching
                          stars in the reference catalog to an image
    @parameter minnbrstars : Minimum number of stars that should match
                          so that a transform is considered valid.
    @parameter maxdist : The maximum shift (in pixels) from the reference
                          catalog to the image to consider

    returns True if all images match
    '''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    print "Ok, we have %i images." % len(db)
    print "Reading alignment star catalog..."
    if refcat is None:
        alistars = readpickle(os.path.join(workdir, "alistars.pkl"))
    else:
        alistars = readpickle(refcat)



    print "Automatic reference image selection :"
    sorteddb = sorted(db, key=lambda k: k['nbrsources'])
    refimg = sorteddb[-1]
    refimgname = refimg['rawimgname']
    print refimgname


    # This is needed as we use a matching that can handle several scales (not yet implemented in alipy)
    refscalingfactor = refimg['scalingfactor']

    # Ok, we are done

    starttime = datetime.now()

    maxnbrautostars = 60
    pairstolerance = tolerance

    for image in db:

            print "- " * 40
            print image['i']+1, "/", len(db), ":", image['rawimgname']

            scalingratio = refimg['scalingfactor']/image['scalingfactor']
            print "scalingratio :", scalingratio

            autostars = readsexcatasstars(image['catfilepath'])
            autostars = sortstarlistbyflux(autostars) # crucial for speed !
            nbrallautostars = len(autostars)
            autostars = autostars[:maxnbrautostars]
            print "Keeping %i stars among %i." % (len(autostars), nbrallautostars)

            (flag, foundangle, foundshift) = findtrans(autostars, alistars, scalingratio, tolerance = tolerance,  minnbrstars = minnbrstars,  mindist = maxdist)

            image['identangle'] = foundangle
            image['foundshift'] = foundshift

            if flag < 0:
                    image['okforalignment'] = False
                    image['alicomment'] = "Cannot be aligned."
                    image['nbralistars'] = 0
                    image['identangle'] = 0.0
                    print "I'll have to skip this one ...\n"
                    continue
            else : 
                    image['okforalignment'] = True
            print "Transformation found."


            (comment, pairs) = formpairs(alistars, autostars, foundangle, foundshift, scalingratio, tolerance = pairstolerance)
            if comment != "":
                    print "Comment :", comment
            image['alicomment'] = comment

            # write the input file for the iraf geomap task
            # "xref yref x y"

            print "Writing geomap input file with %i stars." % len(pairs)
            image['nbralistars'] = len(pairs)

            image['geomapinfile'] = image['basename'] + ".geomapin"	# the filename
            image['geomapinpath'] = os.path.join(workdir, image['geomapinfile'])

            writeforgeomap(image['geomapinpath'], pairs)

            print "Done."
    print "- " * 40



    donetime = datetime.now()
    timetaken = nicetimediff(donetime - starttime)

    print "Here we are. It took me %s" % timetaken


    print "\nSummary :"
    print " -   number of stars that will be used for the alignment."
    print " -   a (rough) rotation angle [deg] that I've found and used for the identification."
    print " -   an eventual comment that should allow you to refine your alignment stars.\n"

    for image in db:
            print "| %4i | %8.3f | %s | %s " % (image['nbralistars'], image['identangle'], image['rawimgname'], image['alicomment'])




    # We "update" the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))

    return reduce(lambda x,y: x and y, [image['okforalignment'] for image
                                        in db])
	

                          
###############################################################

def createTransforms(workdir):
    '''Calls the IRAF routine geomap to construct the pixel by pixel
    tranform for each image'''


    print "\nOK, now let's pyraf !"

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    db = [image for image in db if image['okforalignment'] == True]
    if len(db) == 0:
        print "Cannot continue, no matching images!"
        return
    print "Ok, we have %i images ready for alignment." % len(db)

    # Get the ref image (similar to previous script)

    sorteddb = sorted(db, key=lambda k: k['nbrsources'])
    refimg = sorteddb[-1]
    refimgname = refimg['rawimgname']


    # We get the size of the reference image:
    refimgpath = [image['rawimgpath'] for image in db if image['rawimgname'] == refimgname][0]

    iraf.images()
    iraf.imutil()
    iraf.immatch()

    iraf.unlearn(iraf.imutil)
    imheadblabla = iraf.imutil.imhead(images = image['rawimgpath'], longhea = True, Stdout=1)
    strdim = imheadblabla[0].split("[")[1].split("]")[0].split(",")
    xmax = int(strdim[0])
    ymax = int(strdim[1])

    print "Size of reference image : (%i, %i)" % (xmax, ymax)

    starttime = datetime.now()

    for image in db:

            print "- " * 40
            print image['i']+1, "/", len(db), ":", image['rawimgname']

            image['geodatabasepath'] = os.path.join(workdir, image['basename'] + ".geodatabase")
            if os.path.isfile(image['geodatabasepath']):
                    os.remove(image['geodatabasepath'])
            image['geosumpath'] = os.path.join(workdir, image['basename'] + ".geosum")
            if os.path.isfile(image['geosumpath']):
                    os.remove(image['geosumpath'])

            iraf.unlearn(iraf.geomap)

            iraf.geomap.fitgeom = "rscale"		# You can change this to : shift, xyscale, rotate, rscale
            iraf.geomap.function = "polynomial"	# Surface type
            iraf.geomap.maxiter = 3			# Maximum number of rejection iterations
            iraf.geomap.reject = 3.0		# Rejection limit in sigma
# units



    # other options you could specify :
    #(xxorder=                    2) Order of x fit in x
    #(xyorder=                    2) Order of x fit in y
    #(xxterms=                 half) X fit cross terms type
    #(yxorder=                    2) Order of y fit in x
    #(yyorder=                    2) Order of y fit in y
    #(yxterms=                 half) Y fit cross terms type
    #(calctyp=                 real) Computation type


            iraf.geomap.transfo = "broccoli"	# keep it
            iraf.geomap.interac = "no"		# keep it
            iraf.geomap.verbose = "yes"		# keep it
            iraf.geomap.results = image['geosumpath'] # The optional results summary files

            mapblabla = iraf.geomap(input = image['geomapinpath'], database = image['geodatabasepath'], xmin = 1, xmax = xmax, ymin = 1, ymax = ymax, Stdout=1)

            for line in mapblabla:
                    #print line
                    if "X and Y scale:" in line:
                            mapscale = line.split()[4:6]
                    if "Xin and Yin fit rms:" in line:
                            maprmss = line.split()[-2:]
                    if "X and Y axis rotation:" in line:
                            mapangles = line.split()[-4:-2]
                    if "X and Y shift:" in line:
                            mapshifts = line.split()[-4:-2]

            geomaprms = math.sqrt(float(maprmss[0])*float(maprmss[0]) + float(maprmss[1])*float(maprmss[1]))
            geomapangle = float(mapangles[0])
            geomapscale = float(mapscale[0])

            if mapscale[0] != mapscale[1]:
                    raise mterror("Error reading geomap scale")

            print "Scale :", geomapscale
            print "Angle :", geomapangle
            print "RMS   :", geomaprms

            image['geomaprms'] = geomaprms
            image['geomapangle'] = geomapangle
            image['geomapscale'] = geomapscale


    print "- " * 40

    endtime = datetime.now()
    timetaken = nicetimediff(endtime - starttime)

    print "IRAF geomap done in %s." % timetaken


    print "\nSummary :"
    print " -   geomap RMS"
    print " -   geomap angle [deg]"
    print " -   geomap scale"
    print " -   number of alignment stars.\n"

    for image in db:
            print "| %8.3f | %8.3f | %8.3f | %3i | %s " % (image['geomaprms'], image['geomapangle'], image['geomapscale'], image['nbralistars'], image['rawimgname'])

    print "\nCongrats ! If this worked, you are nearly done."

    # We "update" the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))


############################################################################

def transformImages(workdir):
    '''Tranforms images, based on geomap output from createTransforms, and
    outputs aligned images'''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    db = [image for image in db if image['okforalignment'] == True]
    if len(db) == 0:
        print "Cannot continue, no matching images!"
        return

    print "Ok, we have %i images ready to be transformed." % len(db)

    starttime = datetime.now()

    iraf.images()
    iraf.immatch()

    for image in db:

            print "- " * 40
            print image['i']+1, "/", len(db), ":", image['rawimgname']

            image['aliimgname'] = image['basename'] + "_ali.fits"
            image['aliimgpath'] = os.path.join(workdir, image['aliimgname'])
            image['maskname'] = '%s_mask.fits' % image['basename']

            if os.path.isfile(image['aliimgpath']):
                    os.remove(image['aliimgpath'])


            iraf.unlearn(iraf.gregister)
            iraf.gregister.geometry = "geometric"	# linear, distortion, geometric
            iraf.gregister.interpo = "spline3"	# linear, spline3
            iraf.gregister.boundary = "constant"	# padding with zero
            iraf.gregister.constant = 0.0
            iraf.gregister.fluxconserve = "yes"

            #regblabla = iraf.gregister(input = imgtorotate, output = aliimg, database = databasename, transform = "broccoli", xmin = 1, xmax = dimx, ymin = 1, ymax = dimy, Stdout=1)
            regblabla = iraf.gregister(input = image["rawimgpath"], output =
            image["aliimgpath"], database = image["geodatabasepath"],
            transform = "broccoli", Stdout=1)

            #Prep for coadd: associate a mask file, yet to be made
            hdu = pyfits.open(image['aliimgpath'])[0]
            hdu.header.update('BPM', image['maskname'])
            hdu.writeto(image['aliimgpath'], clobber=True)

            #print regblabla

    print "- " * 40
    endtime = datetime.now()
    timetaken = nicetimediff(endtime - starttime)


    print "Dear user, I'm done with the alignment. I did it in %s." % timetaken

    # We "update" the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))

#############################################################

def transformMasks(workdir, maskext='.mask'):
    '''Apply the image transforms to mask files for each image. If a mask
    file doesn't exist, a mask image of all True values will be used to
    determine the spatial extend of the transformed image.

    A mask should have a 0 where the image is good, and 1 where the image
    is bad. A weightwatcher flag map may be used without problem, but with
    loss of information.
    
    @parameter maskext : If an image file is called image.fits, then its
                 mask file will be called image.maskext.fits
    '''

    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    db = [image for image in db if image['okforalignment'] == True]
    if len(db) == 0:
        print "Cannot continue, no matching images!"
        return

    print "Ok, we have %i images ready for mask transformation." % len(db)

    starttime = datetime.now()

    iraf.images()
    iraf.immatch()

    for image in db:


            print "- " * 40
            print image['i']+1, "/", len(db), ":", image['rawimgname']

            imgdir, imgname = os.path.split(image['rawimgpath'])
            base, ext = os.path.splitext(imgname)
            maskfile = '%s/%s%s%s' % (imgdir, base, maskext, ext)

            if not os.path.exists(maskfile):
                #inputmask='%s/area.fits' % workdir
                inputmask= workdir+'area.fits'

                if not os.path.exists(inputmask):
                    print "1"
                    rawimage = pyfits.open(image['rawimgpath'])[0].data
                    area = np.ones_like(rawimage)
                    hdu = pyfits.PrimaryHDU(area)
                    hdu.writeto(inputmask)

            else:
                print "2"
                #flip the convention so that we can transform the pixels
                mask = pyfits.open(maskfile)[0].data
                mask[mask > 0] = -1
                mask[mask == 0] = 1
                mask[mask == -1] = 0
                mask = mask.astype(np.int16)
                hdu = pyfits.PrimaryHDU(mask)
                
                inputmask = '%s/%s%s%s' % (workdir, base, maskext, ext)
                
                hdu.writeto(inputmask)
                    

            print "3"
            image['maskname'] = image['basename'] + "_mask.fits"
            image['maskpath'] = os.path.join(workdir, image['maskname'])

            if os.path.isfile(image['maskpath']):
                    os.remove(image['maskpath'])


            iraf.unlearn(iraf.gregister)
            iraf.gregister.geometry = "geometric"	# linear, distortion, geometric
            iraf.gregister.interpo = "spline3"	# linear, spline3
            iraf.gregister.boundary = "constant"	# padding with zero
            iraf.gregister.constant = 0.0
            iraf.gregister.fluxconserve = "yes"

            print inputmask
            regblabla = iraf.gregister(input = inputmask, 
                                       output = image["maskpath"], 
                                       database = image["geodatabasepath"], 
                                       transform = "broccoli", Stdout=1)

            print "4"
            #switch the mask back to regular convention
            mask = pyfits.open(image['maskpath'])[0].data
            mask[mask <= 0.95] = -1
            mask[mask > 0.95] = 0
            mask[mask == -1] = 1

            mask = mask.astype('int16')
            hdu = pyfits.PrimaryHDU(mask)
            hdu.writeto(image['maskpath'], clobber=True)
            os.remove(inputmask)

    print "- " * 40
    endtime = datetime.now()
    timetaken = nicetimediff(endtime - starttime)


    print "Dear user, I'm done with mask transformation. I did it in %s." % timetaken

    # We "update" the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))


##############################################################################


def coaddImages(workdir, coaddfile = None, method='median'):
    '''Combine aligned images, masking out bad pixels using transformed
    masks

    @parameter coaddfile path to output coadded file to be produced by
    this method. None defaults to the workdir.
    
    @parameter method Method for combining images. See IRAF's imcombine
    for additional options
    '''
    
    print "Reading database..."
    db = readpickle(os.path.join(workdir, "db.pkl"))
    db = [image for image in db if image['okforalignment'] == True]
    if len(db) == 0:
        print "Cannot continue, no matching images!"
        return

    print "Ok, we have %i images to combine." % len(db)

    curdir = os.getcwd()
    os.chdir(workdir)

    if coaddfile is None:
        coaddfile = 'coadd.fits'

    starttime = datetime.now()

    iraf.images()
    iraf.immatch()

    inputimages = ','.join([image['aliimgname'] for image in db])
    iraf.flprcache()
    iraf.imcombine(input = inputimages, output = coaddfile,
                   combine=method, masktype='goodvalue')

    print "- " * 40
    endtime = datetime.now()
    timetaken = nicetimediff(endtime - starttime)

    os.chdir(curdir)

    print "Dear user, I'm done with the alignment. I did it in %s." % timetaken

    # We "update" the db :
    writepickle(db, os.path.join(workdir, "db.pkl"))


