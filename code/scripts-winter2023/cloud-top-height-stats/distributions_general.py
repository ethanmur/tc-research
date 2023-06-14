# import...
import numpy as np
import os
import sys
import pandas as pd
from geopy import distance
import metpy.calc as mpcalc
import xarray as xr
import math
import matplotlib.pyplot as plt
os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops


# Goal: create a generalized function that finds cloud heights and compiles info for all composite groups (intensity, shear, etc).
# 		do all the sorting later: the information should all be in this dataframe!
# 		can also save locally! easier access later
# 		code is partially based on "code/eye cloud paper figures/Figure 4 cloud dists vs shear quadrant.ipynb"
# Inputs: None! Just run this script to create the nice, organized data file
# 		 No need for a tc='all' input: just do the sorting after making this dataframe!
# Return: a dataframe with height / distribution information for every pass! 
# Notes * -> a list, category = WH, intensifying, DSR depending on input flags! Structure:
# flight | pass | UTC time * | x dists * | y dists * | cloud heights * | vertical bins | normalized dists | defined eyewalls or not | intensity | intensification | tc category | shear strength | shear dir
# case 0
# case 1
# ... 
def find_dists():
	# use a helper fn to get year and file names, along with intensity, etc. metadata
	tc='all'
	yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc=tc)
	metadata = eyewall_metadata.all_metadata()

	# the dataframe below contains simple xy distances and is locally saved. 
	# Recreate this dataset using "code/tests/2023-06-07 find cartesian and shear xy distances.ipynb"
	df_dists = pd.read_pickle("/Users/etmu9498/research/data/simple_distances.pkl")
