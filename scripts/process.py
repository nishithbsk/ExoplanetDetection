# Team Name: DAN
# Members: Arnav Mariwala, Daniel Faniel, Nishith Khandwala

# Calibration Products Rejected
# Biases: 1, 4, 5, 7 (The pixel distribution were not good gaussian curves)
# Darks: None (The pixel counts were similar)
# Flats: None (Their curves were all very similar and good and
#              there were no planes/ satellites running through them.)

#Here are the libraries we need. Note that wherever you see np, that
#stands for Numpy. 'import' loads an external library.

import pyfits
import numpy as np
import sys,os


#Python is an interpreted programming language, so we have to put all of our functions BEFORE
#the main body of the code!


#This function statistically combines bias images into a masterbias image
def AverageBias(biasfiles):
    ''' 
    AverageBias produces a master bias image from a list of individual
    bias exposures.
           biasfiles  - list of file names
           returns a 2D numpy array
    '''

    #opens each bias image file and stores the 2d images in a list
    biasdata= [pyfits.open(i)[0].data for i in open(biasfiles)] 
    
    biascube=np.array(biasdata) #This is a 1-D array of 2-D numpy arrays,
#each element corresponds to a different bias image.

    medianimage=np.median(biascube,axis=0) #The axis=0 tells numpy which
#axis to do the median over. Without this, it would do the median over all
#dimensions of the array and the output would be a single number. With
#axis=0, our output is a 2-D array the same size as an image. 

    meanimage=np.mean(biascube,axis=0)
    sigmaimage=np.std(biascube,axis=0)

    ### !!! TODO FINISH THIS FUNCTION !!!
    masterbias=medianimage

    return masterbias #This is the end of the function 




#This function does the combining of dark currents
def AverageDark(darkfiles,masterbias):

    darkdata=[pyfits.open(i)[0].data for i in open(darkfiles)]
    
    #What is happening here? What is being put into this list?
    # Answer: This list contains the exposure times for each
    # image in the files list, darkfiles.
    darkexpo=[pyfits.open(i)[0].header['exposure'] for i in
              open(darkfiles)]

    darkcube=np.array(darkdata) #This is an array of 2-D images
    darkexpocube=np.array(darkexpo) #This is an array of scalars.
    darklist=[]
    for (image,time) in zip(darkcube,darkexpocube): #The zip command here
        #loops over two lists simultaneously. The iterator of darkcube is
        #image, the iterator of darkexpocube is time

        ### !!! TODO FINISH THIS FUNCTION !!!

        cleandark=(image-masterbias)
        normdark=cleandark/time #(Dark Image- Bias Image) / Time, gives
#image/s which we need, since dark current for any observation depends on
#the exposure time.
        darklist.append(normdark)
    cleandark2=np.array(darklist) 
    masterdark=np.median(cleandark2, axis=0)
    return masterdark


#This function creates a combined flat field image
def AverageFlat(flatfiles,masterbias,masterdark):
    flatdata=[pyfits.open(i)[0].data for i in open(flatfiles)]
    flatexpo=[pyfits.open(i)[0].header['exposure'] for i in
              open(flatfiles)]
    flatcube=np.array(flatdata)
    flatexpocube=np.array(flatexpo)
    cleanflatlist=[]
    for (image,time) in zip(flatcube,flatexpocube):
        cleanflat=image-masterbias-time*masterdark
        normflat=cleanflat/np.median(cleanflat,axis=None)
        cleanflatlist.append(normflat)
    cleancube=np.array(cleanflatlist)
    masterflat=np.median(cleancube,axis=0)
    return masterflat

#This function creates the processed science image after combined bias,
#dark, and flat images have been created.  
def ScienceExposure(rawscidata,masterbias,masterdark,masterflat):
    rawimage=rawscidata.data
    expotime=rawscidata.header['exposure']
    scienceimage=(rawimage-expotime*masterdark-masterbias)/masterflat
    return scienceimage

#Each of these is an argument that needs to be on the calling of the
#script. If the calling does not have 5 additional arguments you will run into errors!

sciencefiles=sys.argv[1]   #First argument is the name of the science observation
flatfilelist=sys.argv[2]  #Second argument is a text file that lists the
#names of all of the flat field images
darkfilelist=sys.argv[3]  #List of dark current image file names
biasfilelist=sys.argv[4] #Fourth is a list of bias image file names
basename=sys.argv[5]    #All of the output files will start with the
#string value of basename. 

finalbias=AverageBias(biasfilelist) 
#Find function above

print 'Bias Image Created'
print finalbias.shape #Change the axis=0 to axis=1 and axis=None and see
#how the shape of the array changes...


finaldark=AverageDark(darkfilelist, finalbias)#What else is needed to get a proper dark current image?)
#Find function above

print 'Dark Current Image Created'
print finaldark.shape #Compare this to the bias image shape. Must these
                      #be the same? What happens if they are not?
                      # Answer: Yes, they must be the same. Otherwise,
                      # numpy will not be able to operate (subtract)
                      # on the images.


finalflat=AverageFlat(flatfilelist,finalbias,finaldark)
#Find function above

for i, sciencefile in enumerate(open(sciencefiles)):
    rawdata=pyfits.open(sciencefile)[0]

    finalimage=ScienceExposure(rawdata,finalbias,finaldark,finalflat)
    #Find function above


    sciheader=rawdata.header #This grabs the header object from the FITS
                          #object rawdata

    newscience=basename+'_'+str(i+1)+'_CleanScience.fits'  #Appending filenames onto the base


    sciencehdu=pyfits.PrimaryHDU(finalimage,header=sciheader)  #This converts
    #a numpy array into a FITS object with a data block (finalimage) and a
    #header (sciheader)

    sciencehdu.writeto(newscience, clobber=True) #This writes the fits object
                          #to the file name newscience, which is defined
                          #above The clobber means to overwrite the file if it
                          #already exists.
                      
newbias=basename+'_Bias.fits'
newdark=basename+'_Master_Dark.fits'
newflat=basename+'_Master_Flat.fits'


biashdu=pyfits.PrimaryHDU(finalbias) #Further writing of files
biashdu.writeto(newbias, clobber=True)

darkhdu=pyfits.PrimaryHDU(finaldark)
darkhdu.writeto(newdark, clobber=True)

#Given the names and conventions, how do we write our final flat field
#image to the disk?
flathdu=pyfits.PrimaryHDU(finalflat)
flathdu.writeto(newflat, clobber=True)