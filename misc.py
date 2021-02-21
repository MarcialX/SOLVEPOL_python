# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Miscellaneous
# Miscellaneous tools
#
# TODO list
# 1. [14/02/2021]. Add Sigma clipping as IDL example
#
# Marcial Becerril, @ 14 February 2021
# Latest Revision: 14 Feb 2021, 01:35 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clip

from scipy.ndimage import median_filter

#from tqdm import tqdm

def sky_level(img, skymode, skysig):

	size = np.shape(img)
	nrow = size[0]
	ncol = size[1]

	npts = nrow*ncol

	maxsky = np.max(2*npts/(nrow-1), 10000)

	# Mantain the same data type as the input image
	istep = npts/maxsky + 1
	nskyvec = maxsky + 200
	skyvec = np.zeros(nskyvec)
	nstep = nrow/istep

	jj = 0
	index0 = istep*np.arange(nstep)
	if nstep > 1:
		i0 = np.max((nrow-1-np.max(index0)-istep)/2, 0)
		index0 = index0 + i0


	for i in range(ncol):
		index = index0 + (i%istep)
		row = img[:,i]

		ng = nrow
		imax = np.max(value_locate(index, [ng-1])[0], 0)
		ix = np.min(index[0:imax], (ng-1))
		skyvec[jj] = row[ix]

		jj = jj + imax + 1 
		if jj >= nskyvec:
			skyvec = [skyvec, skyvec[0:200]*0]

	skyvec = skyvec[0:jj-1]



def value_locate(vector, value):

	vector_diff = np.diff(vector)
	mask = vector_diff > 0
	n = np.sum(mask)
	N = len(vector)

	inc = False
	if n == N-1:
		inc = True
	elif n == 0:
		inc = False
	else:
		return

	result = np.zeros_like(value)
	for i in range(len(value)):
		for j in range(N):
			if inc:
				if value[i] < vector[0]:
					result[i] = -1
					break
				elif (value[i]>=vector[j]) and (value[i]<vector[j+1]):
					result[i] = j
					break 
				elif value[i]>=vector[N-1]:
					result[i] = N-1
					break
			else:
				if vector[0] <= value[i]:
					result[i] = -1
					break
				elif (value[i]>=vector[j+1]) and (value[i]<vector[j]):
					result[i] = j
					break 
				elif value[i]<vector[N-1]:
					result[i] = N-1
					break				

	return result


def sigma_clipping(img, box_width=3, n_sigma=2.5, iterate=5, monitor=True):

	# Make sure width is odd
	box_width = 2*(int(box_width/2))+1
	print("Box width needs to be odd. Box width = ", box_width)

	if box_width < 3:
		print("Box width needs to be greater or equal to 3")
		return

	if iterate < 1 or iterate > 20:
		print("Iterations between 1 - 20")
		return 

	bw = box_width**2

	mean = (median_filter(img, size=box_width)*bw - img)/(bw - 1)

	if n_sigma <= 0:
		print("Sigma must be greater than 0")
		return 

	imdev = (img - mean)**2
	fact = (n_sigma**2)/(bw - 2)

	imvar = fact*(median_filter(imdev, size=box_width)*bw - imdev)

	wok = np.where(imdev < imvar)
	nok = wok[0].size

	npix = img.size
	nchange = npix - nok

	if monitor:
		print('{:.2f}% of pixels replaced, n_sigma = {:.1f}'.format(nchange*100./npix, n_sigma))

	# Bad pixels are equal to the full pixels
	if nok == npix:
		print("All the pixels can't be replaced")
		return

	if nok > 0:
		mean[wok] = img[wok]

	# Solve for n iterations
	if n_sigma >= 2:
		iterate -= 1
		if iterate == 0:
			return mean 
		return sigma_clipping(mean, box_width=box_width, n_sigma=n_sigma, iterate=iterate, monitor=monitor)

	return mean