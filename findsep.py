# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Find separation
# Find separation in pixels between ordinary and extraordinary sources
#
# Input:
# 	filename: Path of the ordinary/extraordinary image
#
# Outputs:
#	x_off, y_off: Separation in pixels between the images
#
# Marcial Becerril based on Edgar Andre SOLVEPOL repo
# https://github.com/edgarandre/SOLVEPOL/blob/master/findsep.pro, 
# @ 21 August 2021
# Latest Revision: 21 Aug 2021, 22:00 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #


import os

import numpy as np

from astropy.io import fits
from scipy.signal import correlate2d
from scipy import ndimage

from misc import *

from matplotlib.pyplot import *
ion()


def xcorr2d(x, window='hamming'):
	"""
		Cross-correlation
		Parameters
		----------
		x : 2d-np.array
			Image to get the cross-correlation
		window: string
			Windowing method
		----------
	"""
	# Image size
	x_size, y_size = np.shape(x)

	# If image has equals dimensions
	if x_size == y_size and not window is None:
		# Windowing method
		if window == 'blackman':
			w = np.abs(np.blackman(x_size))
		elif window == 'hamming':
			w = np.abs(np.hamming(x_size))
		elif window == 'hanning':
			w = np.abs(np.hanning(x_size))
		else:
			print('Windowing method not defined')
			return
		# Get 2D window
		window = np.sqrt(np.outer(w, w))

		# Windowing the image 
		x = x * window

	# Get the FFT of the image
	x = np.fft.fft2(x)

	# Get the conjugate in order to calculate the correlation
	xs = np.conj(x)

	# Compute the conjugate times the FFT and then the inverse FFT
	# Cross-correlation according to the Convolution theorem
	corr = np.fft.fftshift(np.fft.ifft2(xs*x))

	return np.abs(corr)


def find_separation(img, margin=0.0, mask_size=(10, 10), window='hamming'):
	"""
		Find separation in pixels between ordinary/extraordinary images
		Parameters
		----------
		img : fits image
			Polarimetric image
		margin : float
			Margin images
		mask_size : tuple
			Mask size for the center
		----------
	"""
	# Get image data
	img_data = img[0].data
	# Get size image
	img_size = np.shape(img_data)

	# Get dimensions
	x_size = img_size[0]
	y_size = img_size[1]

	# Define border
	x_border = int(margin*x_size)
	y_border = int(margin*y_size)

	# Get cross-correlation
	xcorr = xcorr2d(img_data, window=window)

	# Mask the center region
	if x_size%2 == 1:
		x_center = int(y_size/2) + 1
	else:
		x_center = int(y_size/2)

	if y_size%2 == 1:
		y_center = int(y_size/2) + 1
	else:
		y_center = int(y_size/2)

	mask_zeros = np.zeros((2*mask_size[0], 2*mask_size[1]))
	xcorr[x_center-mask_size[0]:x_center+mask_size[0], y_center-mask_size[1]:y_center+mask_size[1]] = mask_zeros

	# Get maximum correlation position
	max_corr = np.sort(np.ravel(xcorr))[-1]
	max_pos = np.where(xcorr == max_corr)

	if np.shape(max_pos)[1] > 1:
		max_pos = max_pos[0]
	else:
		max_pos = [max_pos[0][0], max_pos[1][0]]

	x_off = max_pos[0]-x_center
	y_off = max_pos[1]-y_center

	# Print offset in pixels
	print ("=======================================")
	print("X separation: ", x_off, " pixels")
	print("Y separation: ", y_off, " pixels")
	print ("=======================================")

	return x_off, y_off


# Defining initial parameters
# =================================================
margin = 0.0
mask_size = (10, 10)

# Read FITS file
# =================================================
filename = "./HD126593_L9_0_corrected.fits"
img = fits.open(filename)

xpos, ypos = find_separation(img, margin=margin, mask_size=mask_size)
