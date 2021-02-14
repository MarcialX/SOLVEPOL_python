# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Reduction
# Data-Reduction tools
#
# Input:
# 	objectfile: text file listing image names
#	flats: text file listing flat names
#	bias: text file listing bias name
#	subdir: subdirectory where the reduced images will be located
#
# Outputs:
#	Several FITS files:
#		*bias_zero.fits
#		*flat_bias_combine.fits
#		*bias_flat_objnames.fits
#
# TODO list
# 1. [07/02/2021]. Save corrected images
#
# Marcial Becerril, @ 24 January 2021
# Latest Revision: 14 Feb 2021, 01:28 GMT-6
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

from tqdm import tqdm

from matplotlib.pyplot import *

ion()


def get_fits_from_list(file_list_path):
	"""
		Get Fits Images from a List
		Parameters
		----------
		bias_file : string
			Text file path of the file which contains the fits data
		----------
	"""
	# Get list from the txt file
	fits_list = {}
	with open(file_list_path) as files:
		for file in files:
			file = file.split('\n')[0]
			# Get name
			name = file.split('/')[-1]
			name = name.split('.')[0]
			fits_list[name] = {}
			# File name 
			fits_list[name]['path'] = file
			# Read fits file
			fits_list[name]['fits'] = fits.open(file)


	return fits_list


def combine_images(int_img, size_x, size_y, method):
	"""
		Combine images
		Parameters
		----------
		int_img : data cube
			Data cube: l x m x n
			l : number of images
			m : x-size
			n : y-size
		size_x : int
			x length
		size_y : int
			y length
		method : string
			Combination method
		----------
	"""
	comb = np.zeros((size_x, size_y))
	
	# Combining the images from the cube
	# =================================================
	for i in tqdm(range(size_x), desc='Combining images'):
		for j in range(size_y):
			if method == 'avg':
				comb[i,j] = np.mean(int_img[:,i,j]) 
			elif method == 'med':
				comb[i,j] = np.median(int_img[:,i,j]) 

	return comb


def overscan_polynom(int_img, x_ovr_scn, size):
	"""
		Fit overscan with a polynom
		Parameters
		----------
		int_img : array
			image to extract the polynom
		x_ovr_scan : list
			x_ovr_scn[0] : True if x-axis, False if y-axis
			x_ovr_scn[1:2] : min and max limits of the scan 
		size : list
			size[0] : x-size
			size[1] : y-size
		----------
	"""
	xoverscan = x_ovr_scn[0]
	ovr_scn_1 = x_ovr_scn[1]
	ovr_scn_2 = x_ovr_scn[2]

	overscan_flag = False

	# Get size
	size_x = size[0]
	size_y = size[1]

	axis_x = np.arange(size_x)
	axis_y = np.arange(size_y)

	pols_img = np.zeros((size_x, size_y))

	# Correct flat images by overscan
	# =============================================
	if xoverscan:
		# Over X-AXIS
		# =============================================
		img_overscan = np.zeros(size_x)

		for m in range(size_x):
			img_range = int_img[m,ovr_scn_1:ovr_scn_2]
			if len(img_range) > 0:
				overscan_flag = True
				img_overscan[m] = np.mean(img_range)

		if overscan_flag:
			# Fitting the overscan region of the bias
			fx = np.polyfit(axis_x, img_overscan, 2)
			p = np.poly1d(fx)

			# Array with polynom
			for n in range(size_y):
				pols_img[:,n] = p(axis_x)

	else:
		# Over Y-AXIS
		# =============================================
		img_overscan = np.zeros(size_y)

		for n in range(size_y):
			img_range = int_img[ovr_scn_1:ovr_scn_2,n]
			if len(img_range) > 0:
				overscan_flag = True
				img_overscan[n] = np.mean(img_range)

		if overscan_flag:
			# Fitting the overscan region of the bias
			fy = np.polyfit(axis_y, img_overscan, 2)
			p = np.poly1d(fy)	

			# Array with polynom
			for m in range(size_x):
				pols_img[m,:] = p(axis_y)

	return pols_img


def write_fits(name, header, data):
	"""
		Get bias corrected image
		Parameters
		----------
		name : string
			Name, with the full path, of the file
		header : header
			Fits header
		data : data
			Fits data 
		----------
	"""
	os.system('rm -rf '+name+'.fits')
	hdu = fits.PrimaryHDU(list(data), header=header)
	hdu.writeto(name+'.fits')


def get_corrected_bias(file_list_path, x_ovr_scn, ref='bias_0001', comb_method='avg', sig_clip=[5, 5], save=True):
	"""
		Get bias corrected image
		Parameters
		----------
		file_list_path : string
			Path of file with the list of bias fits images paths
		x_ovr_scan : list
			x_ovr_scn[0] : True if x-axis, False if y-axis
			x_ovr_scn[1:2] : min and max limits of the scan 
		ref : string
			Reference image. 'bias_0001' by default 
		comb_method : string
			Combination method: 'avg' or 'med'
		sig_clip : list
			sig_clip[0] : sigma threshold
			sig_clip[1] : max number of iterations
		----------
	"""
	# Read the images
	fits_list = get_fits_from_list(file_list_path) 

	# Reference image
	ref_img = fits_list[ref]['fits'][0]
	size_x = ref_img.header['NAXIS1']
	size_y = ref_img.header['NAXIS2']
	size = [size_x, size_y]
	noImgs = len(fits_list.keys())

	int_img = np.zeros((noImgs, size_x, size_y))
	pol_img = np.zeros((size_x, size_y))

	# Data cube
	# =================================================
	for i, img in enumerate(fits_list.keys()):
		int_img[i,:,:] = fits_list[img]['fits'][0].data

	# Combining bias images
	zero_com = combine_images(int_img, size_x, size_y, comb_method)

	# Averaging over all columns of the overscan region
	# =================================================
	# Get polynom
	pol_img = overscan_polynom(zero_com, x_ovr_scn, size)

	# Sigma clipping
	# =================================================
	sigma = sig_clip[0]
	maxiters = sig_clip[1] 
	zero_com = sigma_clip(zero_com, sigma=sigma, maxiters=maxiters)

	# Get the bias zero
	bias_zero = zero_com - pol_img

	# Saving FITS file
	# =================================================
	if save:
		name = "bias_corrected"
		print ("Saving FITS bias file as: "+name)
		print ("=======================================")
		write_fits(name, ref_img.header, bias_zero)

	print("Done")

	return bias_zero, zero_com


def get_corrected_flats(file_list_path, x_ovr_scn, bias_zero, ref='flat_0001', comb_method='avg', sig_clip=[5, 5], save=True):
	"""
		Get flats corrected image
		Parameters
		----------
		file_list_path : string
			Path of file with the list of flats fits images paths
		x_ovr_scan : list
			x_ovr_scn[0] : True if x-axis, False if y-axis
			x_ovr_scn[1:2] : min and max limits of the scan 
		bias_zero : 2d-array
			Bias image corrected
		ref : string
			Reference image. 'flat_0001' by default 
		comb_method : string
			Combination method: 'avg' or 'med'
		sig_clip : list
			sig_clip[0] : sigma threshold
			sig_clip[1] : max number of iterations
		----------
	"""
	# Read the images
	fits_list = get_fits_from_list(file_list_path)

	# Reference image
	ref_img = fits_list[ref]['fits'][0]
	size_x = ref_img.header['NAXIS1']
	size_y = ref_img.header['NAXIS2']
	size = [size_x, size_y]
	noImgs = len(fits_list.keys())

	int_img = np.zeros((size_x, size_y))
	median_flat = np.zeros((size_x, size_y))
	flats_bias = np.zeros((size_x, size_y))
	
	pols_img = np.zeros((noImgs, size_x, size_y))
	flat_no_bias = np.zeros((noImgs, size_x, size_y))

	# Reading data cube
	# =================================================
	i = 0
	for img in tqdm(fits_list.keys(), desc='Correcting overscan'):

		int_img = fits_list[img]['fits'][0].data

		# Averaging over all columns of the overscan region
		# =================================================
		# Get polynom
		pols_img[i,:,:] = overscan_polynom(int_img, x_ovr_scn, size)

		# Get Flat-Bias
		# =============================================
		flat_no_bias[i,:,:] = int_img - bias_zero - pols_img[i,:,:]
		i += 1

	# Combining bias images
	flat_comb = combine_images(flat_no_bias, size_x, size_y, comb_method)

	# Sigma clipping
	# =================================================
	sigma = sig_clip[0]
	maxiters = sig_clip[1] 
	#flat_comb = sigma_clip(flat_comb, sigma=sigma, maxiters=maxiters)

	# Saving FITS file
	# =================================================
	if save:
		name = "flats_corrected"
		print ("Saving FITS flats file as: "+name)
		print ("=======================================")
		write_fits(name, ref_img.header, flat_comb)

	return flat_comb


def correct_images(file_list_path, bias_corrected, flats_corrected, x_ovr_scn, ref="HD126593_L0_1"):
	"""
		Correct images
		Parameters
		----------
		file_list_path : string
			Path of file with the list of objects fits images
		bias_corrected : 2d-array
			Bias corrected image
		flats_corrected : 2d-array
			Flats corrected image
		x_ovr_scan : list
			x_ovr_scn[0] : True if x-axis, False if y-axis
			x_ovr_scn[1:2] : min and max limits of the scan
		ref : string
			Reference image. 'HD126593_L0_1' by default 
		----------
	"""
	# Read the images
	fits_list = get_fits_from_list(file_list_path)

	# Reference image
	ref_img = fits_list[ref]['fits'][0]
	size_x = ref_img.header['NAXIS1']
	size_y = ref_img.header['NAXIS2']
	size = [size_x, size_y]
	noImgs = len(fits_list.keys())

	int_img = np.zeros((size_x, size_y))
	poly_array = np.zeros((noImgs, size_x, size_y))
	
	# BIAS CORRECTION
	# =================================================

	print ("Correcting by bias")
	print ("=======================================")

	object_no_bias = np.zeros((noImgs, size_x, size_y))

	# Reading data cube of each object
	# =================================================
	i = 0
	for img in tqdm(fits_list.keys(), desc='Correcting by bias'):
		
		int_img = fits_list[img]['fits'][0].data

		# Averaging over all columns of the overscan region
		# =================================================
		# Get polynom
		poly_array[i,:,:] = overscan_polynom(int_img, x_ovr_scn, size)

		# Get Flat-Bias
		# =============================================
		object_no_bias[i,:,:] = int_img - bias_corrected - poly_array[i,:,:]
		i += 1

	print ("Done")

	# FLAT CORRECTION
	# =================================================

	print ("Correcting by flats")
	print ("=======================================")

	object_no_flats = np.zeros((noImgs, size_x, size_y))

	# Flat correction
	# =================================================
	medValFlat = np.median(flats_corrected)
	for i in tqdm(range(noImgs), desc='Correcting by flats'):
				
		object_no_flats[i,:,:] = (object_no_bias[i,:,:])*(medValFlat/ flats_corrected)

	print ("Done")

	return object_no_flats



subdir = 'solvepol'


# CCD data: 
# NOTA. Hay varios modelos con valores distintos.
# Agregar mini base de datos o JSON file

# Andor iKon-L936: 2048 X 2048
# No overscan
xoverscan = False
overscan1 = 0
overscan2 = -1

# ===================== B I A S =====================
# Sigma clipping parameters
sig_clip = [5, 5]

overscan = [xoverscan, overscan1, overscan2]
bias_zero, bias_comb = get_corrected_bias('./HD126593/bias.list', overscan, sig_clip=sig_clip)

# ==================== F L A T S ====================
flat_comb = get_corrected_flats('./HD126593/flats.list', [xoverscan, 0, 0], bias_zero, sig_clip=[5, 5])

# ======== C O R R E C T I N G   I M A G E S ========
objects = correct_images('./HD126593/stdstar.list', bias_zero, flat_comb, [xoverscan, 0, 0])


# Objectos reales
# fits_list = get_fits_from_list('./HD126593/stdstar.list')
# img = list(fits_list.keys())[0]
# obj_obs = fits_list[img]['fits'][0].data


# NOTAS
# 1. flat_comb no debería restarsele el bias zero, sería mejor quitarle 
# el zero_com.
# 2. Como esta funcionando el overscan (Posibles errores en el IDL)
