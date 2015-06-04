'''
photometry.py

This is a stub file provided for lab 1.3 part I. 

COPY THIS FILE TO A WORK DIRECTORY OF YOUR CHOOSING FIRST!!!!!!

Fill in the missing blanks in the following functions.
'''

###############################################

import numpy as np
import fitmodel
import pyfits, ldac
import sys
from sextractor import *
import matplotlib.pyplot as plt
import os
import string

################################################
# 1.0528
def aperFlux(source_radius = 10, background_width = 10, 
             gain = 5.752, nimages = 1, errflag=True, image=None, x=None, y=None,
             exptime = 0., darkexptime = 0.,
             masterbias = None, masterdark = None):
    ''' aperFlux - Measures the flux (in adc counts), gaussian flux
             uncertainty, and the probability that the background is
             described by a constant.
         @parameter image - A 2D numpy array containing the image, from
             hdu.data
         @parameter x - The x position of the object, in pixels
         @parameter y - The y position of the object, in pixels
         @parameter source_radius - The radius of a circular region
             centered on the object representing the source region
         @parameter background_width - Width of an annulus, with
             source_radius < R <= source_radius + background_width,
             representing the background region
         @parameter gain - The gain of the detector
         @parameter nimages - The number of images that were coadded to
             form this image
         @parameter errflag- This bool determines whether or not errors
             should be determined for the flux. Only used starting with
             Lab 1.5, once we figure out the calculation of uncertainties

         @returns flux, fluxerr, background_prob
         '''

    ysize, xsize = image.shape

    #create 2D arrays with the X,Y coordinates of every pixel
    X,Y = np.meshgrid(np.arange(xsize), np.arange(ysize))

    #calculate the distance of every pixel from the object
    #note: the x-1 accounts for the offset between sextractor and array indicies
    dR = np.sqrt((X - x - 1)**2 + (Y - y - 1)**2)

    #2D boolean arrays selecting the background and source regions
    inSource = dR <= source_radius
    inBackground = np.logical_and(dR > source_radius, 
                                  dR <= (background_width + source_radius))
    inRegion = dR <= source_radius + background_width

    #counting pixels
    nsourcepix = len(image[inSource])
    nbackgroundpix = len(image[inBackground])
    npixels = len(image[inRegion])

    ### TODO: FINISH THESE LINES
    # calculate the flux of the source, the flux uncertainty, and the
    # probability that the background is described by a constant
    #
    # Feel free to add additional intermediate steps as needed. We will
    # need to calculate
    # the flux with Lab 1.3, but the uncertainties on the flux will wait
    # until Lab 1.5.  Set the fluxerr and background_prob to zero until
    # then. 

    source_image = sum(image[inSource]) * gain
    background_image = sum(image[inBackground]) * gain

    if nbackgroundpix:
        avg_background = float(background_image) / nbackgroundpix
    else:
        avg_background = 0

    flux = float(source_image) - (nsourcepix * avg_background)  

    # Readout
    bias_image = (sum(masterbias[inSource]) * gain) #######
    avg_bias = bias_image / nsourcepix

    var_bias, background_chi2 = 0, 0
    iter_start = 0 - source_radius - background_width
    iter_end = source_radius + background_width
    for x_ep in xrange(iter_start, iter_end):
        x_modified = x + x_ep
        if x_modified < 0 or x_modified >= xsize: continue
        for y_ep in xrange(iter_start, iter_end):
            y_modified = y + y_ep
            if y_modified < 0 or y_modified >= ysize: continue

            if inRegion[y_modified][x_modified]:
                var_bias += np.square((masterbias[y_modified][x_modified] * gain) - avg_bias)
            if inBackground[y_modified][x_modified]:
                background_chi2 += np.square((image[y_modified][x_modified] * gain) - 
                                   avg_background) / (image[y_modified][x_modified] * gain)


    readout = np.sqrt(nimages) * np.sqrt(var_bias / (npixels - 1))
    # DarkCurr
    dark_image = (sum(masterdark[inSource]) * gain * nimages) / darkexptime
    dark_image_avg = dark_image / npixels


    ####

    f_obj = flux * nimages / exptime
    f_sky = background_image * nimages / exptime
    sig = f_obj * exptime
    noise = np.sqrt(f_sky*exptime*nsourcepix + f_obj*exptime +
                    dark_image_avg*exptime*nsourcepix + nsourcepix*(readout**2))

    sig_to_noise = sig / noise
    N = len(image[inBackground]) - 1
    background_prob = fitmodel.chisq_exceeds_prob(background_chi2, N)
    # fluxerr, background_prob = 0, 0

    # if (errflag):
    #     fluxerr = np.sqrt(abs(flux))
    #     background_chi2 = 0.
    #     for count in image[inBackground]:
    #         background_chi2 += ((count - avg_background)**2) / count
    #     dof = nbackgroundpix - 1
    #     background_prob = fitmodel.chisq_exceeds_prob(background_chi2, dof)

    fluxerr = noise
    return flux, fluxerr, background_prob

#############################################

def aperMag(flux, fluxerr, errflag=False):
    ''' aperMag - converts a measured flux and fluxerr into a magnitude
    and magnitude error. Assumes gaussian error propagation.
         @parameter flux  the flux of the object (may be a 1D array)
         @parameter fluxerr the gaussian flux uncertainty of the object
                          (may be a 1D array)
         @returns mag, magerr  - the magnitude and gaussian uncertainty in
                                    the magnitude
    '''

    ### TODO: FINISH THESE LINES
    # calculate the magnitude and the gaussian magnitude uncertainty.
    #
    # Feel free to add additional intermediate steps as needed. Raw
    # magnitudes will be calculated starting in Lab 1.3, and errors will
    # be propagated in Lab 1.5. Set the magerr term to zero until then. 
    if flux < 0: flux = 1
    mag = -2.5 * np.log10(flux) 
    magerr = fluxerr * np.sqrt(((-2.5 / np.log(10))/flux)**2)

    return mag, magerr
    
    
###############################################

def createPhotCat(image, detectcat, masterbias, masterdark):
    ''' createPhotCat - Given an image and a detection catalog, this
    function will calculate the flux at every x,y position in the
    detection catalog and return the results as a new catalog
         @parameter image - A pyfits HDU containing a coadded image
         @parameter detectcat - An ldac catalog containing an xcol and a
                                ycol
         @parameter otherparams - Any additional function arguments passed
    to this function with key=val pairs will be passed to aperFlux.

         @returns an LDAC catalog with a flux, fluxerr, mag, magerr, and
         background_prob column
         '''

##### TODO:  Either 1.) Use this function and write your own to combine
##### the photcats from an R, G, & B exposure into one catalog.
##### or
##### 2.) Rewrite this function to do photometry on RGB at put it into one catalog

    nimages = image.header['NCOMBINE']
    exptime = image.header['exposure']
    darkexptime = image.header['exposure']

    #define some arrays to hold our calculations
    fluxs = np.zeros(len(detectcat))
    fluxerrs = np.zeros(len(detectcat))
    backprobs = np.zeros(len(detectcat))
    mags = np.zeros(len(detectcat))
    magerrs = np.zeros(len(detectcat))

    #loop over the objects to measure
    for i, x,y in zip(np.arange(len(detectcat)),
                      detectcat['X_IMAGE'], 
                      detectcat['Y_IMAGE']):

        #calculate the flux; propagate parameters to the aperFlux function
        fluxs[i], fluxerrs[i], backprobs[i] = \
            aperFlux(image=image.data, x=x, y=y, nimages = nimages, exptime=exptime, 
                     darkexptime=darkexptime, masterbias=masterbias.data,
                     masterdark=masterdark.data)

        mags[i], magerrs[i] = aperMag(fluxs[i], fluxerrs[i])

    # return fluxs, fluxerrs, backprobs, mags, magerrs

    #create columns for a new catalog
    newcols=[]
    newcols.extend([pyfits.Column(name = 'flux',
                                  format = 'E',
                                  array = fluxs),
                    pyfits.Column(name = 'fluxerr',
                                  format = 'E',
                                  array = fluxerrs),
                    pyfits.Column(name = 'backprob',
                                  format = 'E',
                                  array = backprobs),
                    pyfits.Column(name = 'mag',
                                  format = 'E',
                                  array = mags),
                    pyfits.Column(name = 'magerr',
                                  format = 'E',
                                  array = magerrs)])

    #convert the columns into an LDAC catalog, copying the existing
    #columns from the detection cat
    photcat = ldac.LDACCat(pyfits.new_table(detectcat.hdu.columns + \
                                                pyfits.ColDefs(newcols), 
                                            header=detectcat.hdu.header))
    return photcat

        
##############
masterbias_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/normal_data/process_output/normal_Master_Bias.fits'
masterdark_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/normal_data/process_output/normal_Master_Dark.fits'

masterdark = pyfits.open(masterdark_file_path)[0]
masterbias = pyfits.open(masterbias_file_path)[0]

base_name = sys.argv[1]

folder = '/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/normal_data/reduction_output/coadd_images/'
all_coadd_images = os.listdir(folder)
for i, image in enumerate(all_coadd_images):
	full_path = folder + image
	obs = pyfits.open(full_path)[0]
	output_base_name = base_name + '_' + str(i)
	ref_catalog_file = output_base_name + '_refdetect.cat'
	ref_file = str("/afs/ir.stanford.edu/class/physics100/workdir/DAN_Project/data/normal_data/photometry_output/ref_files/" + ref_catalog_file)
	sextractor(full_path, ref_file)
	ref_catalog = ldac.openObjectFile(ref_file)

	# New Lab 1.5
	catalog = createPhotCat(obs, ref_catalog, masterbias, masterdark)

	flux = catalog['flux']
	fluxerr = catalog['fluxerr']
	mag = catalog['mag']
	magerr = catalog['magerr']
	backprob = catalog['backprob']


	catalog.saveas('/usr/class/physics100/workdir/DAN_Project/data/normal_data/photometry_output/normal_files/' + output_base_name + '.cat', clobber=True)

