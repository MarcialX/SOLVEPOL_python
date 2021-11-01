# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-Runphot
# Aperture photometry of the paired sources
#
# Input:
#
# Outputs:
#
# Marcial Becerril based on Edgar Andre SOLVEPOL repo
# https://github.com/edgarandre/SOLVEPOL/blob/master/dither.pro,
# @ 17 October 2021
# Latest Revision: 17 Oct 2021, 21:36 GMT-6
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


def read_pair_file(filename):
	"""
		Read pair source file
		Parameters
		----------
		filename : str
			Input file containing pairing file
		----------
	"""
	# Open file
	with open(filename, 'r') as file:
		data = file.readlines()

	pairs = {}
	obj_cnt = 1
	flag_pair = False

	for line in data:
		# Get pair name
		if 'PAIR ' in line:
			pairnum = line.split('PAIR ')[1]
			pair_name = 'P'+pairnum[:-1]
			if pairnum[:-1].isnumeric():
				pairs[pair_name] = {}
				obj_cnt = 1
				flag_pair = True

		elif flag_pair:
			line_fields = line.split('\t')
			line = []
			for i in range(len(line_fields)):
				lf = line_fields[i].replace(' ', '')
				if not lf == '':
					line.append(lf)

			# Check if line has the enough columns
			if len(line) >= 6 and line[0].isnumeric():
				obj_name = 'OBJ'+str(obj_cnt)
				pairs[pair_name][obj_name] = [int(line[0]), float(line[1]), float(line[2]),\
											  float(line[3]), float(line[4]), float(line[5])]
				obj_cnt += 1

	return pairs


# Parameters
# Photons per Analog Digital Units
PHPADU = 1.
# Photometry aperture radii
APR = 10.
# Inner and outer radii to be used for the sky annulus
SKYRAD = (5, 10)
# Two element vector giving the minimum and maximum value
# of a good pixel
BADPIXEL = (25, 50)


# Read fits image
filename = './HD126593_L0_0_corrected.fits'
header, data = read_fits_file(filename)

# Read pair sources file
pairname = './sources.pair'
pair_objs = read_pair_file(pairname)
