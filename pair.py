# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Pairing
# Search for the pair objects: ordinary-extraordinary
#
# Input:
#
# Outputs:
#
# Marcial Becerril based on Edgar Andre SOLVEPOL repo
# https://github.com/edgarandre/SOLVEPOL/blob/master/dither.pro, 
# @ 11 September 2021
# Latest Revision: 11 Sep 2021, 00:03 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

import os
import time

import numpy as np

from astropy.io import fits

from tqdm import tqdm

from misc import *
from fits_tools import read_fits_file

from matplotlib.pyplot import *
ion()


# CONSTANTS
MAX_PAIRS = 32767 

def pairing(filename, shift, delta=(2,2), deltamag=1, **kwargs):
	"""
		Find the pair objects given the shift distance between
		ordinary and extraordinary sources
		Parameters
		----------
		filename : str
			Input file containing results of finding
		shift : tuple
			Shift of two polarization components along the x-axis in pixels
		delta : tuple
			The +/- x/y-axis search range
		deltamag : float
			The magnitude difference limit for finding pairs
		**kwargs:
			verbose : bool
				Show internal information
		----------
	"""

	# Key arguments
	# ----------------------------------------------
	# Write file
	write_file = kwargs.pop('write_file', True)
	# Name file
	pair_filename = kwargs.pop('pair_filename', 'sources.pair')
	# Verbose key
	verbose = kwargs.pop('verbose', False)
	# ----------------------------------------------

	# Extract shift x and y
	shiftx = shift[0]
	shifty = shift[1]

	# Delta x/y axis range
	deltax = delta[0]
	deltay = delta[1]

	# Read founded sources file
	source_file = np.loadtxt(filename, skiprows=1)

	# Read header
	with open(filename, 'r') as f:
	    header_list = f.readline()[:-1]

	# Convert to list
	header_list = header_list.split('\t')

	# Parse data as dictionary
	src_data = {}
	for i, param in enumerate(header_list):
		src_data[param] = source_file[:,i]

	# Get # rows and columns
	nlines, cols = source_file.shape

	# Get flux in magnitudes
	mag = 25 - 2.5*np.log10(src_data['FLUX'])

	# Initialise pairs counter
	pairnum = 0
	nopairs = 0
	nmagcut = 0

	if write_file:
		# Create pair file
		pair_file = open(pair_filename, 'w')
		# Get system time
		systime = time.asctime(time.localtime(time.time()))
		pair_file.write('\nProgram: PAIR '+systime+'\n')
		pair_file.write('\nPairing criteria:\nshiftx = '+str(shiftx)+' shifty = '+str(shifty)+\
						'\ndeltax = '+str(deltax)+' deltay = '+str(deltay)+'\n')
		pair_file.write('\nNo. of sources rejected by SPACIAL criteria')
		pair_file.write('\nNo. of sources rejected by DELTAMAG criteria\n\n')
		pair_file.write('  \tSTAR\tX\tY\tFLUX\tSHARP\tROUND\n')

	if verbose:
		print ('  \tSTAR\tX\tY\tFLUX\tSHARP\tROUND\n')

	pairs_objs = {}

	for i in range(nlines):
		alpha = src_data['X'][i] + shiftx
		beta = src_data['Y'][i] + shifty 

		valid_points = np.where((src_data['X'] > alpha-deltax) & (src_data['X'] < alpha+deltax) & \
								(src_data['Y'] > beta-deltay) & (src_data['Y'] < beta+deltay))[0]

		if len(valid_points) > 0:
			valid_mag_points = np.where(np.abs(mag[valid_points]-mag[i]) < deltamag)[0]
			if len(valid_mag_points) > 0:
				k = valid_points[valid_mag_points]
				best_point = np.argmin((src_data['X'][k]-alpha)**2+(src_data['Y'][k]-beta)**2)
				j = k[best_point]
				pairnum += 1

				pair_name = 'P'+str(pairnum)
				pairs_objs[pair_name] = {}
				pairs_objs[pair_name]['OBJ1'] = [src_data['STAR'][i], src_data['X'][i], src_data['Y'][i],\
							src_data['FLUX'][i], src_data['SHARP'][i], src_data['ROUND'][i]]
				pairs_objs[pair_name]['OBJ2'] = [src_data['STAR'][j], src_data['X'][j], src_data['Y'][j],\
							src_data['FLUX'][j], src_data['SHARP'][j], src_data['ROUND'][j]]

				if write_file:
					pair_file.write('PAIR '+str(pairnum)+'\n')
					pair_file.write('  \t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'%(src_data['STAR'][i], src_data['X'][i], src_data['Y'][i],\
							src_data['FLUX'][i], src_data['SHARP'][i], src_data['ROUND'][i]))
					pair_file.write('  \t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'%(src_data['STAR'][j], src_data['X'][j], src_data['Y'][j],\
							src_data['FLUX'][j], src_data['SHARP'][j], src_data['ROUND'][j]))
				if verbose:
					print('PAIR '+str(pairnum)+'\n')
					print('  \t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'%(src_data['STAR'][i], src_data['X'][i], src_data['Y'][i],\
							src_data['FLUX'][i], src_data['SHARP'][i], src_data['ROUND'][i]))
					print('  \t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'%(src_data['STAR'][j], src_data['X'][j], src_data['Y'][j],\
							src_data['FLUX'][j], src_data['SHARP'][j], src_data['ROUND'][j]))			
			else:
				nmagcut += 1

	nopairs = nlines - 2*pairnum - nmagcut

	pair_file.close()

	# Display the pairing results
	print ('No. of pairs = ', pairnum)
	print ('No. of non-pairs = ', nopairs)

	if pairnum > MAX_PAIRS:
		print ('The number of pairs is greater than the maximum allowed ('+str(MAX_PAIRS)+')')
		return

	if pairnum == 0:
		print ('There are no sources (pairs) left.')
		print ('Verify x-y shift, and/or try with lower flux/sigma')
		return

	if write_file:
		# Edit file updating pairs founded
		# ----------------------------------------------
		with open(pair_filename, 'r') as file:
			data = file.readlines()

		# Modifiy the rejected results
		data[7] = 'No. of sources rejected by SPACIAL criteria: '+str(nopairs)+'\n'
		data[8] = 'No. of sources rejected by DELTAMAG criteria: '+str(nmagcut)+'\n\n'

		# Update file
		with open(pair_filename, 'w') as file:
			file.writelines(data)
		# ----------------------------------------------

	return pairs_objs


# EXAMPLE
# ----------------------------------------------
# Define pairing parameters
filename_pairs = './sources.find'
# Shift obtained from findsep.py
shift = (-57, 57)
# Get pairs
pairs = pairing(filename_pairs, shift, delta=(3,3), deltamag=1)

# Display image
filename_img = './HD126593_L0_0_corrected.fits'
header, data = read_fits_file(filename_img)

imshow(np.log10(data), origin='lower', vmin=1.85, vmax=2.3)
xticks([],[])
yticks([],[])

for pair in pairs.keys():
	# Coordinates object 1
	obj1x = pairs[pair]['OBJ1'][1]
	obj1y = pairs[pair]['OBJ1'][2] 
	# Coordinates object 2
	obj2x = pairs[pair]['OBJ2'][1]
	obj2y = pairs[pair]['OBJ2'][2] 

	plot(obj1x, obj1y, 'rx')
	plot(obj2x, obj2y, 'co')