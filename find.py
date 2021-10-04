# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Find stars
# Find positive brightness perturbations (i.e stars) in an image
# 	Also returns centroids and shape parameters (roundness & sharpness).
# 	Adapted from 1991 version of DAOPHOT, but does not allow for bad pixels
# 	and uses a slightly different centroid algorithm.
#   Modified in March 2008 to use marginal Gaussian fits to find centroids  
#
# Input:
# 	Image - 2 dimensional image array (integer or real) for which one
#	wishes to identify the stars present
#
# Outputs:
#
# Marcial Becerril based on IDL find.pro
# 
# @ 11 September 2021
# Latest Revision: 11 Sep 2021, 00:03 GMT-6
# https://idlastro.gsfc.nasa.gov/ftp/pro/idlphot/find.pro
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import os

import numpy as np

from scipy.signal import convolve2d
from astropy.io import fits

from misc import *
from fits_tools import *

from tqdm import tqdm

from matplotlib.pyplot import *
ion()


def find_objs(data, *args, **kwargs):
	"""
		Find brightness perturbations in an image, e.g. stars
		Parameters
		----------
		data : np.array-2d
			Image to find brigntess objects
		*args:
			hmin : float
				Threshold intensity for a point source - should generally 
				be 3 or 4 sigma above background RMS
			fwhm : float
				Full-width at half maximum (in pixels)
			sharplim : 2-items tuple
				2 element vector giving low and high cutoff for the
				sharpness statistic
			roundlim : 2-item tuple
				2 element vector giving low and high cutoff for the
				roundness statistic
		**kwargs:
			verbose : bool
				Show internal information
		----------
	"""

	# Arguments
	# ----------------------------------------------
	# Get total number of parameters
	npar = len(args)
	# Initial arguments
	init_args = [10., 1., (0.2, 1.0), (-1.0, 1.0)]

	for i, arg in enumerate(args):
		init_args[i] = arg

	# Threshold
	hmin = init_args[0]
	# Full-width at half maximum
	fwhm = init_args[1]
	# Sharpness limits
	sharplim = init_args[2]
	# Roundness limits
	roundlim = init_args[3]
	# ----------------------------------------------

	# Key arguments
	# ----------------------------------------------
	# Verbose key
	verbose = kwargs.pop('verbose', True)
	# ----------------------------------------------

	# Initial parameters
	# -----------------------------------------
	# Maximum size of convolution box in pixels
	maxbox = 13


	# Check if data loaded is 2-dimensional
	img_size = np.shape(data)
	if len(img_size) != 2:
		print('ERROR - Image array must be 2 dimensional')
		return

	# Asign dimensions
	size_x = img_size[0]
	size_y = img_size[1]
	if verbose:
		print('Input image size is '+str(size_x)+', '+str(size_y))

	# Check FWHM level
	if fwhm < 0.5:
		print('WARNING - Supplied FWHM must be at least 0.5 pixels')
		print('FWHM set in 0.5 pixels')
		fwhm = 0.5

	# Get radius, equals to 1.5 sigma
	radius = max(0.637*fwhm, 2.001)
	rad_sq = radius**2
	nhalf = min(int(radius), int((maxbox-1)/2))
	# Number of pixels inside of convolution box
	nbox = 2*nhalf + 1
	# Index of central pix
	middle = nhalf

	lastro = size_x - nhalf
	lastcl = size_y - nhalf
	sigsq = (fwhm/2.35482)**2
	
	# Gaussian convolution kernel container
	g = np.zeros((nbox, nbox))

	dd = np.arange(nbox-1, dtype=int) + 0.5 - middle
	dd2 = dd**2

	row2 = (np.arange(nbox, dtype=float)-nhalf)**2

	for i in range(nhalf+1):
		temp = row2 + i**2
		g[nhalf-i] = temp
		g[nhalf+i] = temp

	# Mask identifies valid pixels in convolution box
	mask = np.where(g < rad_sq, 1, 0)
	good = np.where(g < rad_sq)
	# Number of valid pixels
	pixels = len(good[0])

	# Compute quantities for centroid computations that can be used for all stars
	g = np.exp(-0.5*g/sigsq)

	# Getting the weights
	xwt = np.zeros((nbox, nbox))
	wt = nhalf - np.abs(np.arange(nbox, dtype=float)-nhalf) + 1 
	p = np.sum(wt)
	for i in range(nbox):
		xwt[i] = wt

	ywt = xwt.T

	sgx = np.sum(g*xwt, axis=1)
	sgy = np.sum(g*ywt, axis=0)

	sumgx = np.sum(wt*sgy)
	sumgy = np.sum(wt*sgx)

	sumsqy = np.sum(wt*sgy**2)
	sumsqx = np.sum(wt*sgx**2)

	vec = nhalf - np.arange(nbox, dtype=float)

	dgdx = sgy*vec
	dgdy = sgx*vec

	sdgdxs = np.sum(wt*dgdx**2)
	sdgdx = np.sum(wt*dgdx)

	sdgdys = np.sum(wt*dgdy**2)
	sdgdy = np.sum(wt*dgdy)

	sgdgdx = np.sum(wt*sgy*dgdx)
	sgdgdy = np.sum(wt*sgx*dgdy)

	# Convolution kernel now in c
	c = g*mask 
	sumc = np.sum(c)
	sumcsq = np.sum(c**2) - sumc**2/pixels
	sumc = sumc/pixels
	c[good] = (c[good]-sumc)/sumcsq

	c1 = np.exp(-0.5*row2/sigsq)
	sumc1 = np.sum(c1)/nbox
	sumc1sq = np.sum(c1**2) - sumc1
	c1 = (c1-sumc1)/sumc1sq

	if verbose:
		print('RELATIVE ERROR computed from FWHM '+str(np.sqrt(np.sum(c[good]**2))))

	# CONVOLUTION
	# ======================================
	print('Beginning convolution of image')

	h = convolve2d(data, c, mode='same')

	minh = np.min(h)
	h[:nhalf], h[size_x-nhalf:size_x] = minh, minh 
	h[:,:nhalf], h[:,size_y-nhalf:size_y] = minh, minh

	# Convert to 1d array
	h = h.flatten()

	print('Convolution of image finished')
	# ======================================

	# Excluding central pixel
	mask[middle, middle] = 0
	# Updating information
	pixels = pixels - 1
	good = np.where(mask)

	# x-y coordinate of valid pixels relative to the center
	xx = np.mod(good[0], nbox) - middle
	yy = np.mod(good[1], nbox) - middle

	offset = yy*size_x + xx

	# Valid pixel greater than hmin
	print(hmin)
	index = np.array(np.where(h >= hmin)).flatten()
	nfound = len(index)

	if nfound == 0:
		print('ERROR - No maxima exceed input threshold of '+str(hmin))
		return 

	for i in range(pixels):
		
		stars = np.where(h[index] >= h[index+offset[i]])
		nfound = len(stars[0])
		if nfound == 0:
			print('ERROR - No maxima exceed input threshold of '+str(hmin))
			return	
		index = index[stars]

	ix = np.mod(index, size_x)				# X index of local maxima
	iy = np.array(index/size_y, dtype=int)	# Y index of local maxima

	ngood = len(index)

	nstar = 0		# NSTAR counts all stars meeting selection criteria
	badround = 0
	badsharp = 0
	badcntrd = 0

	# Create output X and Y arrays?
	if npar >= 2:
		x = np.zeros(ngood)
		y = np.copy(x)

	# Create output flux, sharpness arrays?
	if npar >= 4:
		flux = np.copy(x)
		sharp = np.copy(x)
		roundness = np.copy(x)

	# Create output file?
	if verbose:
		print('Create output file...')

	# Create results table
	print('STAR\tX\tY\tFLUX\tSHARP\tROUND')

	# Loop over star positions.
	# Compute statistics

	for i in range(ngood):
		temp = data[iy[i]-nhalf:iy[i]+nhalf+1, ix[i]-nhalf:ix[i]+nhalf+1]
		d = h[index[i]]

		# Compute sharpness statistic
		sharp1 = (temp[middle, middle]-(np.sum(mask*temp))/pixels)/d
		if sharp1 < sharplim[0] or sharp1 > sharplim[1]:
			badsharp += 1
			continue

		# Compute roundness
		dx = np.sum(np.sum(temp, axis=0)*c1)
		dy = np.sum(np.sum(temp, axis=1)*c1)
		if (dx <= 0) or (dy <= 0):
			badround += 1
			continue

		around = 2*(dx-dy)/(dx+dy)
		if around < roundlim[0] or around > roundlim[1]:
			badround += 1
			continue

		# Centroid computation
		# ==============================================
		# Find X centroid
		sd = np.sum(temp*ywt, axis=0)

		sumgd = np.sum(wt*sgy*sd)
		sumd = np.sum(wt*sd)
		sddgdx = np.sum(wt*sd*dgdx)

		# Height of the best-fitting marginal Gaussian. If this is not
		# positive then the centroid does not make sense
		hx = (sumgd - sumgx*sumd/p) / (sumsqy - sumgx**2/p)

		if hx < 0:
			badcntrd += 1
			#continue

		skylvl = (sumd - hx*sumgx)/p 
		dx = (sgdgdx - (sddgdx-sdgdx*(hx*sumgx + skylvl*p)))/(hx*sdgdxs/sigsq)
		if np.abs(dx) >= nhalf:
			badcntrd += 1
			continue

		# X centroid in original array
		xcen = ix[i] + dx

		# Find Y centroid
		sd = np.sum(temp*xwt, axis=1)

		sumgd = np.sum(wt*sgx*sd)
		sumd = np.sum(wt*sd)
		sddgdy = np.sum(wt*sd*dgdy)

		hy = (sumgd-sumgy*sumd/p)/(sumsqx - sumgy**2/p)

		if hy < 0:
			badcntrd += 1
			continue

		skylvl = (sumd - hy*sumgy)/p 
		dy = (sgdgdy - (sddgdy-sdgdy*(hy*sumgy + skylvl*p)))/(hy*sdgdys/sigsq)
		if np.abs(dy) >= nhalf:
			badcntrd += 1
			continue

		ycen = iy[i] + dy 

		# This star has met all selection criteria. Print out and save results
		if verbose:
			print ('%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f'%(nstar, xcen, ycen, d, sharp1, around))

		if npar >= 2:
			x[nstar] = xcen
			y[nstar] = ycen

		if npar >= 4:
			flux[nstar] = d
			sharp[nstar] = sharp1
			roundness[nstar] = around

		nstar += 1

	# NSTAR is now the index of the last star found
	nstar -= 1

	if verbose:
		print('No. of sources rejected by sharpness criteria: ', badsharp)
		print('No. of sources rejected by roundness criteria: ', badround)
		print('No. of sources rejected by centroid criteria: ', badcntrd)

	if nstar > 0:
		if npar >= 2:
			x = x[:nstar+1]
			y = y[:nstar+1]

		if npar >= 4:
			flux = flux[:nstar+1]
			sharp = sharp[:nstar+1]
			roundness = roundness[:nstar+1]

		if npar >= 4:
			return x, y, flux, sharp, roundness
		elif npar >= 2:
			return x, y

	else:
		print('No source founded')
		return



# Read FITS file
filename = './HD126593_L0_0_corrected.fits'
header, data = read_fits_file(filename)

x, y = find_objs(data, 30, 10)