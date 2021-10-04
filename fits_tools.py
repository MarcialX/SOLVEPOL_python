# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------- #
# SOLVEPOL-FITS Tools
# Some FITS tools
#
# Marcial Becerril based on Edgar Andre SOLVEPOL repo
# https://github.com/edgarandre/SOLVEPOL/blob/master/dither.pro, 
# @ 11 September 2021
# Latest Revision: 11 Sep 2021, 00:20 GMT-6
#
# For all kind of problems, requests of enhancements and bug reports, please
# write to me at:
#
# mbecerrilt92@gmail.com
# mbecerrilt@inaoep.mx
#
# --------------------------------------------------------------------------------- #

from astropy.io import fits


def read_fits_file(file_path):
	"""
		Load FITS image
		Parameters
		----------
		file_path : string
			Fits file path
		----------
	"""
	# Read the FITS image
	fits_img = fits.open(file_path)
	# Get header
	fits_header = fits_img[0].header
	# Get data
	fits_data = fits_img[0].data

	return fits_header, fits_data
