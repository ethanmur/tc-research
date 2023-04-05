import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import scipy
import sys
import warnings
import pandas as pd
from scipy.signal import find_peaks

os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns
# import helpful_stats


# like other flight level datasets, calculate the mean cloud height vs rmw!
# use the crl data instead of flight level data, though
def calc_cloud_heights_rmws( tc='all', binwidth=.1, maxbin=10):
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"
    yearlist, filelist = helper_fns.get_crl_datasets( tc)
    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of crl files to be added: " + str( filecount))


    # initialize a new dataframe to save binned data
    # do this outside the main loop to store all values in the same dataframe
    df_new_bins = pd.DataFrame( )

    # make an array representing the bins!
    bins = np.arange(0, maxbin, step=binwidth)
    # create an array of midpoints (right between two given binwidths)
    midpoints = []
    for i in range( len( bins) - 1):
        midpoints.append( (bins[ i] + bins[ i+1]) / 2)
    midpointarray = np.array( midpoints)

    # add midpoints and bins to the data array
    df_new_bins[ 'bins'] = bins[ 0: len( bins)-1] # drop last value for correct array size
    df_new_bins[ 'midpoints'] = midpoints
    nancount = 0

    # initial value, update later!
    total_vars = []

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            # get data
            crl_path = crl_data_root + yearval
            os.chdir( crl_path)
            crl_data = xr.open_dataset( fileval)

            print( crl_path)
            print( fileval)
