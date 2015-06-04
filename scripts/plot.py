import math
import numpy as np
import ldac
import matplotlib.pyplot as plt
import sys
import os
import re 

object_data = {
    'normal': {
    	'path': '/usr/class/physics100/workdir/DAN_Project/data/normal_data/photometry_output/',
        'catalog': '/usr/class/physics100/workdir/DAN_Project/data/normal_data/photometry_output/normal_files/',
        'true_mag': 9.42,
        'err_mag': 0.06,
        'refstar_index': 12,
    },
    'eclipse': {
    	'path': '/usr/class/physics100/workdir/DAN_Project/data/eclipse_data/photometry_output/',
        'catalog': '/usr/class/physics100/workdir/DAN_Project/data/eclipse_data/photometry_output/eclipse_files/',
        'true_mag': 9.42,
        'err_mag': 0.06,
        'refstar_index': 12,
    }
}

# MAX ERROR
ERROR_THRESHOLD = 0.5

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def get_error_bar(object_name, cat, corr):
	y_vals = cat['mag'] + corr['corr']
	y_errs = np.sqrt(np.square(cat['magerr']) + corr['partial_err'])

	return_y, return_y_err = [], []
	for y, y_err in zip(y_vals, y_errs):
		if math.isnan(y) is False and math.isnan(y_err) is False:
			if y_err < ERROR_THRESHOLD:
				return_y.append(y)
				return_y_err.append(y_err)

	return (return_y, return_y_err)
	
def plot(object_name, y, y_errs):
	x = np.arange(len(y))
	plt.errorbar(x, y, yerr=y_errs, fmt='o')
	plt.gca().invert_yaxis()
	plt.xlabel("Time")
	plt.ylabel("Magnitude")
	plt.show()

def corrections(refcatalog, data):
    i = data['refstar_index']
    m = {}
    m['corr'] = data['true_mag'] - refcatalog['mag'][i]
    m['partial_err'] = data['err_mag']**2 + refcatalog['magerr'][i]**2
    return m

object_name = sys.argv[1]
data = object_data[object_name]
all_files = os.listdir(data['catalog'])
all_files = sorted_nicely(all_files)

all_ys, all_y_errs = [], []
for cat_file in all_files:
	cat_file = data['catalog'] + cat_file
	cat = ldac.openObjectFile(cat_file)

	refcat = cat

	corr = corrections(refcat, data)

	ys, y_errs = get_error_bar(object_name, cat, corr)
	for y, y_err in zip(ys, y_errs):
		all_ys.append(y)
		all_y_errs.append(y_err)

plot(object_name, all_ys, all_y_errs)
