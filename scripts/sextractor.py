'''
sextractor.py

This library is a python wrapper around the Linux program Source
Extractor.

Source Extractor finds objects in images and creates catalogs of those
objects.

Any extra parameters to source extractor may be passed as keyword
arguments to functions in this library.
'''

import subprocess, os

#########################################################


sextractor_prog="/afs/ir/class/physics100/src/ldacpipeline-0.12.20/bin/Linux_64/sex_theli"
alipy_sex_configdir = '/afs/ir.stanford.edu/class/physics100/src/alipy_1.3/sextractor_config'

##############################################################

def sextractor(image, catalog_name, detect = None, config = None, 
               callMethod = subprocess.check_call,
               **keywords):
    '''
    Calls the program sextractor from python.
    @param image (str - filename) image to measure photometry from
    @param outputcat (str - filename) path and name for output catalog
    @param detect (str - filename) detection image to use if using dual
               image mode
    @param config (str - filename) configuration file to use, otherwise default
    @param callMethod (function) for debugging purposes -- function that makes the system call

    All keywords are appended to the command line in source extractor format. Values that are lists are 
       converted to comma seperated lists in the system call.
    '''
    # print catalog_name
    image = os.path.abspath(image)
    catalog_name = os.path.abspath(catalog_name)

    configFlag = ''
    if config:
        configFlag = '-c %s' % os.path.abspath(config)

    imageFlag = image
    if detect:
        detect = os.path.abspath(detect)
        imageFlag = '%s,%s' % (detect,image)


    command = '%s %s %s' % (sextractor_prog, configFlag, imageFlag)

    if 'catalog_name' not in keywords:
        keywords['catalog_name'] = catalog_name

    if 'catalog_type' not in keywords:
        keywords['catalog_type'] = 'FITS_LDAC'


    isiterable = lambda obj: hasattr(obj,'__iter__')

    params = sorted(keywords.keys())
    for param in params:
        val = keywords[param]
        if isiterable(val):
            val=','.join(map(str, val))
        elif isinstance(val, type(True)):
            if val:
                val = 'Y'
            else:
                val = 'N'
        else:
            val = str(val)
        command = command + ' -%s %s' % (param.upper(), val)

    curdir = os.getcwd()
    os.chdir(alipy_sex_configdir)

    try:

        callMethod(command.split())

    finally:
        os.chdir(curdir)

    return command
