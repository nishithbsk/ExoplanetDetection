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
import sextractor
import matplotlib.pyplot as plt


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
if sys.argv[2] == "m3":  
    red_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/data_reduction_workdir_m3/red_data_reduction_workdir/coadd.fits'
    green_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/data_reduction_workdir_m3/green_data_reduction_workdir/coadd.fits'
    blue_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/data_reduction_workdir_m3/blue_data_reduction_workdir/coadd.fits'
    ref_catalog = ldac.openObjectFile('/afs/ir.stanford.edu/class/physics100/workdir/DAN/ref_detect_m3.cat')
    red_masterbias = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_m3/red_out_Bias.fits'
    red_masterdark = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_m3/red_out_Master_Dark.fits'
    green_masterbias = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_m3/green_out_Bias.fits'
    green_masterdark = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_m3/green_out_Master_Dark.fits'
    blue_masterbias = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_m3/blue_out_Bias.fits'
    blue_masterdark = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_m3/blue_out_Master_Dark.fits'

if sys.argv[2] == "ngc2158":
    red_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/data_reduction_workdir_NGC2158/red_data_reduction_workdir/coadd.fits'
    green_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/data_reduction_workdir_NGC2158/green_data_reduction_workdir/coadd.fits'
    blue_file_path = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/data_reduction_workdir_NGC2158/blue_data_reduction_workdir/coadd.fits'
    ref_catalog = ldac.openObjectFile('/afs/ir.stanford.edu/class/physics100/workdir/DAN/ref_detect_NGC2158.cat')
    red_masterbias = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_NGC2158/red_out_Bias.fits'
    red_masterdark = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_NGC2158/red_out_Master_Dark.fits'
    green_masterbias = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_NGC2158/green_out_Bias.fits'
    green_masterdark = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_NGC2158/green_out_Master_Dark.fits'
    blue_masterbias = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_NGC2158/blue_out_Bias.fits'
    blue_masterdark = '/afs/ir.stanford.edu/class/physics100/workdir/DAN/output_NGC2158/blue_out_Master_Dark.fits'

red, green, blue = pyfits.open(red_file_path)[0], pyfits.open(green_file_path)[0], pyfits.open(blue_file_path)[0]
red_masterdark, green_masterdark, blue_masterdark = pyfits.open(red_masterdark)[0], pyfits.open(green_masterdark)[0], pyfits.open(blue_masterdark)[0]
red_masterbias, green_masterbias, blue_masterbias = pyfits.open(red_masterbias)[0], pyfits.open(green_masterbias)[0], pyfits.open(blue_masterbias)[0]

ims_rgb = [red, green, blue]
output_base_name = sys.argv[1]

# Our group uses the Red filter.
# Generally, we have observed the red filter to give us the better baseline.
# Stars generate more light in the red wavelength.
# ref_catalog_file = output_base_name + '_detect.cat'
# sextractor.sextractor(red, ref_catalog_file)

#ref_catalog = ldac.openObjectFile('/afs/ir.stanford.edu/class/physics100/workdir/DAN/ref_detect_NGC2158.cat')

# Old Lab 1.5
# always, r, then g, then b
# new_cols = []
# mags_rgb = []
# for i in xrange(len(ims)):
#     fluxs, fluxerrs, backprobs, mags, magerrs = createPhotCat(ims_rgb[i], ref_catalog)
#     mags_rgb.append(mags)
#     name = str(i)
#     new_cols.extend([
#         pyfits.Column(name = name + 'flux', format='E', array=fluxs),
#         pyfits.Column(name = name + 'fluxerr', format='E', array=fluxerrs),
#         pyfits.Column(name = name + 'backprob', format='E', array=backprobs),
#         pyfits.Column(name = name + 'mag', format='E', array=mags),
#         pyfits.Column(name = name + 'magerr', format='E', array=magerrs),
#     ])

# blue_minus_green = mags_rgb[2] - mags_rgb[1]
# green_minus_red = mags_rgb[1] - mags_rgb[0]
# new_cols.extend([
#     pyfits.Column(name = 'b_minus_g', format='E', array=blue_minus_green),
#     pyfits.Column(name = 'g_minus_r', format='E', array=green_minus_red),  
# ])

# photcat = ldac.LDACCat(pyfits.new_table(
#     ref_catalog.hdu.columns + pyfits.ColDefs(new_cols), 
#     header=ref_catalog.hdu.header
# ))

# photcat.saveas(outfile_base_name + '.cat', clobber=True)

# New Lab 1.5
red_catalog = createPhotCat(red, ref_catalog, red_masterbias, red_masterdark)
green_catalog = createPhotCat(green, ref_catalog, green_masterbias, green_masterdark)
blue_catalog = createPhotCat(blue, ref_catalog, blue_masterbias, blue_masterdark)

colors_blue_green = blue_catalog['mag'] - green_catalog['mag']
colors_green_red = green_catalog['mag'] - red_catalog['mag']
blue_green_err = np.sqrt(blue_catalog['magerr']**2 + green_catalog['magerr']**2)
green_red_err = np.sqrt(green_catalog['magerr']**2 + red_catalog['magerr']**2)

fluxR, fluxV, fluxB = red_catalog['flux'], green_catalog['flux'], blue_catalog['flux']
fluxerrR, fluxerrV, fluxerrB = red_catalog['fluxerr'], green_catalog['fluxerr'], blue_catalog['fluxerr']
magR, magV, magB = red_catalog['mag'], green_catalog['mag'], blue_catalog['mag']
magerrR, magerrV, magerrB = red_catalog['magerr'], green_catalog['magerr'], blue_catalog['magerr']
backprobR, backprobV, backprobB = red_catalog['backprob'], green_catalog['backprob'], blue_catalog['backprob']


new_cols = []
new_cols.extend([pyfits.Column(name='fluxR',format='E',array=fluxR),
                 pyfits.Column(name='fluxV',format='E',array=fluxV),
                 pyfits.Column(name='fluxB',format='E',array=fluxB),
                 pyfits.Column(name='B - V',format='E',array=colors_blue_green),
                 pyfits.Column(name='V - R',format='E',array=colors_green_red),
                 pyfits.Column(name='BVerr',format='E',array=blue_green_err),
                 pyfits.Column(name='VRerr',format='E',array=green_red_err),
                 pyfits.Column(name='fluxerrR',format='E',array=fluxerrR),
                 pyfits.Column(name='fluxerrB',format='E',array=fluxerrB),
                 pyfits.Column(name='fluxerrV',format='E',array=fluxerrV),
                 pyfits.Column(name='magR',format='E',array=magR),
                 pyfits.Column(name='magB',format='E',array=magB),
                 pyfits.Column(name='magV',format='E',array=magV),
                 pyfits.Column(name='magerrR',format='E',array=magerrR),
                 pyfits.Column(name='magerrB',format='E',array=magerrB),
                 pyfits.Column(name='magerrV',format='E',array=magerrV),
                 pyfits.Column(name='backprobR',format='E',array=backprobR),
                 pyfits.Column(name='backprobB',format='E',array=backprobB),
                 pyfits.Column(name='backprobV',format='E',array=backprobV)])

combined_catalog = ldac.LDACCat(pyfits.new_table(green_catalog.hdu.columns + pyfits.ColDefs(new_cols), header=green_catalog.hdu.header))

combined_catalog.saveas(output_base_name + '_combined.cat', clobber=True)


# OLD CODE !!!!
# plt.gca().invert_yaxis()
# plt.scatter(combined_catalog['Blue - Green'], combined_catalog['mag'])
# plt.title('Mag Vs. Blue - Green')
# plt.ylabel('mag')
# plt.xlabel('blue-green')
# plt.show()

# plt.gca().invert_yaxis()
# plt.scatter(combined_catalog['Green - Red'], combined_catalog['mag'])
# plt.title('Mag Vs. Green - Red')
# plt.ylabel('mag')
# plt.xlabel('green - red')
# plt.show()

# plt.scatter(combined_catalog['mag'], combined_catalog['Green - Red'])
# plt.title('Green - Red Vs. Mag')
# plt.xlabel('mag')
# plt.ylabel('Green - Red')
# plt.show()
