import math
import numpy as np
import ldac
import matplotlib.pyplot as plt
import sys

object_data = {
    'm3': {
        'catalog': 'files/m3.cat',
        'refcatalog': 'files/m3_ref.cat',
        'true_b_mag': 10.89,
        'err_b_mag': 0.06,
        'true_v_mag': 10.13,
        'err_v_mag': 0.03,
        'refstar_index': 71,
    },
    'ngc2158': {
        'catalog': 'files/ngc2158.cat',
        'refcatalog': 'files/ngc2158_ref.cat',
        'true_b_mag': 10.07,
        'err_b_mag': 0.03,
        'true_v_mag': 9.98,
        'err_v_mag': 0.04,
        'refstar_index': 252,
    }
}

# MAX ERROR
ERROR_THRESHOLD = 1.5

def get_error_bar(object_name, cat, corr, x_filter_1, x_filter_2, y_filter):
	# x_filters for plotting over x axis
	# y filters for plotting over y axis
	
	x1_corr = cat['mag' + x_filter_1] + corr['corr' + x_filter_1]
	x2_corr = cat['mag' + x_filter_2] + corr['corr' + x_filter_2]
	x1_err  = np.sqrt(np.square(cat['magerr' + x_filter_1]) + corr['partial_err' + x_filter_1])
	x2_err  = np.sqrt(np.square(cat['magerr' + x_filter_2]) + corr['partial_err' + x_filter_2])
	x_vals  = x1_corr - x2_corr
	x_errs  = np.sqrt(np.square(x1_err) + np.square(x2_err))

	y_vals  = cat['mag' + y_filter] + corr['corr' + y_filter]
	y_errs  = np.sqrt(np.square(cat['magerr' + y_filter]) + corr['partial_err' + y_filter])

	return_x, return_y = [], []
	return_x_err, return_y_err = [], []
	for x, x_err, y, y_err in zip(x_vals, x_errs, y_vals, y_errs):
		# check for NaN - baaad code
		if math.isnan(x) is False and math.isnan(x_err) is False and math.isnan(y) is False and math.isnan(y_err) is False:
			# cutting off error
			# print x_err, y_err
			if x_err < ERROR_THRESHOLD and y_err < ERROR_THRESHOLD:
				# print "YAY"
				return_x.append(x)
				return_x_err.append(x_err)
				return_y.append(y)
				return_y_err.append(y_err)

	print return_x
	return (return_x, return_y, return_x_err, return_y_err)
	
def plot(object_name, x_filter_1, x_filter_2, y_filter, data):
	x, y, x_errs, y_errs = data
	plt.errorbar(x, y, xerr=x_errs, yerr=y_errs, fmt='o')
	plt.gca().invert_yaxis()
	plt.xlabel(x_filter_1 + ' - ' + x_filter_2)
	plt.ylabel(y_filter)
	plt.show()

def corrections(refcatalog, data):
    i = data['refstar_index']
    m = {}
    m['corrV'] = data['true_v_mag'] - refcatalog['magV'][i]
    m['partial_errV'] = data['err_v_mag']**2 + refcatalog['magerrV'][i]**2
    m['corrB'] = data['true_b_mag'] - refcatalog['magB'][i]
    m['partial_errB'] = data['err_b_mag']**2 + refcatalog['magerrB'][i]**2
    return m

object_name = sys.argv[1]
data = object_data[object_name]
cat = ldac.openObjectFile(data['catalog'])
refcat = ldac.openObjectFile(data['refcatalog'])
corr = corrections(refcat, data)
err_data = get_error_bar(object_name, cat, corr, 'B', 'V', 'V')
plot(object_name, 'B', 'V', 'V', err_data)
