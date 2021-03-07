# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Dither
# Combine images
#
# Explanation:
#
# The procedure find the shift between dithered images
# and  combine them. The output images have the
# size of the first image per position of the waveplate. 
#
# Input:
# 	imgsfile : Name of text file containing the names of
#			   the images to combine in a column.
#	ndither	 : Number of images per waveplate position.
#
# Outputs:
#	One combined image per waveplate position.
#
# TODO list
# 1. [07/03/2021]. Dither function missing
# 2. [07/03/2021]. Order functions
#
# Marcial Becerril based on Edgar Andre SOLVEPOL repo
# https://github.com/edgarandre/SOLVEPOL/blob/master/dither.pro, 
# @ 7 March 2021
# 
# Latest Revision: 7 March 2021, 00:28 GMT-6
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
from astropy.stats import sigma_clip

from reduction import *
from misc import *

from tqdm import tqdm

from matplotlib.pyplot import *
ion()


def img_comb(img_list, shift, nwp, ref=None, *args, **kwargs):
	"""
		Combine Fits Images from a List
		Parameters
		----------
		img_list : string
			File path which contains the fits images paths.
		shift	 : array
			x and y shift
		nwp   : int
			Waveplate position
		*args : additional arguments
		**kwargs : additional keywords (for verbose)
		----------
	"""
	# Check for verbose
	verbose = kwargs.pop('verbose', None)

	# Read the images
	fits_list = get_fits_from_list(img_list) 

	# Reference image
	if ref is None:
		ref = list(fits_list.keys())[0]
	ref_img = fits_list[ref]['fits'][0]
	size_x = ref_img.header['NAXIS1']
	size_y = ref_img.header['NAXIS2']
	size = [size_x, size_y]
	noImgs = len(fits_list.keys())

	int_img = np.zeros((noImgs, size_x, size_y))
	
	# Print data if verbose
	if verbose:
		print("******************************************")
		print("Images: ")
		for f in fits_list.keys():
			print(fits_list[f]['path'])
		print("Number of images: ", noImgs)
		print("Waveplate position: ", nwp)
		print("******************************************")

	# Read the images
	# =================================================
	for i, img in enumerate(fits_list.keys()):
		int_img[i,:,:] = fits_list[img]['fits'][0].data

	# Taking the biggest size images
	# =================================================
	xshift = shift[0]
	yshift = shift[1]

	max_xshift = np.max(np.abs(xshift))
	max_yshift = np.max(np.abs(yshift))

	new_int_img = np.zeros((noImgs, size_x+2*max_xshift, size_y+2*max_yshift))
	img_combine = np.zeros((size_x+2*max_xshift, size_y+2*max_yshift))

	# Shifting image reference to first image
	# =================================================
	for i in tqdm(range(noImgs), desc='Shifting images'):
		for j in range(size_x):
			for k in range(size_y):
				# Shift Assignment
				new_int_img[i, j+max_xshift+xshift[i], k+max_yshift+yshift[i]] = int_img[i,j,k] 

	# Adding images
	# =================================================
	for j in range(size_x+2*max_xshift):
		for k in range(size_y+2*max_yshift):

			img_combine[j, k] = np.mean(new_int_img[:,j,k])



	return new_int_img



path_fits_list = './HD126593/stdstar.list'
n = 4
shift = [[1,2,3],[1,2,3]]

d = img_comb(path_fits_list, shift, n, verbose=True)

