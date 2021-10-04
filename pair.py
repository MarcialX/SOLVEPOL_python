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

import numpy as np

from astropy.io import fits

from misc import *

from tqdm import tqdm

from matplotlib.pyplot import *
ion()


subdir = 'solvepol'

ext = '.cut'

#readcol