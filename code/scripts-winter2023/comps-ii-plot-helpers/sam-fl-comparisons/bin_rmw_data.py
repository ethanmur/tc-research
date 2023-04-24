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

os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns

sys.path.append("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023
import helpful_stats

# use the new rmw axis saved in the processed fl data to bin important quantities!
# return the binned quantities for additional analysis :)
# do this for specific
def bin( tc='all', binwidth=.1, maxbin=10, new_data=False, bintype='rmw'):

    if new_data:
        fl_data_root = "/Users/etmu9498/research/data/in-situ-noaa-trimmed/"
    else:
        fl_data_root = "/Users/etmu9498/research/data/in-situ-noaa-processed/"
    yearlist, filelist = helper_fns_winter2023.get_fl_datasets( tc)
    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))

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

            print(filei)
            print(fileval)

            # get data
            fl_path = fl_data_root + yearval
            os.chdir( fl_path)
            fl_data = xr.open_dataset( fileval)

            # bin all of the available variables!
            # input variables
            #input_vars = [ 'dist', 'rmw', 'float_time', 'WS.d', 'UWZ.d', 'SfmrRainRate.1', 'THETAE.d', 'MR.d', 'TA.d'] # change me!
            #var_names = [ 'dist', 'rmw', 'time', 'wind_speed', 'w', 'Rain Rate', 'Theta E', 'Mixing Ratio', 'temp']  # change me!
            
            # input_vars = list( fl_data.keys())
            
            input_vars = ['center_dist', 'rmw', 'WS.d', 'TA.d', 'MR.d', 'UWZ.d', 'THETAE.d', 'SfmrRainRate.1'] 

            ################
            ## the list below holds about 10 empty lists, one for each of the input_vars.
            ## for each variable, there are len( midpoints) empty lists to hold each sorted data point
            ## structure
            ################

            print(input_vars)

            # new 3/17/23 code
            # initial case: make a bunch of empty lists. If not true, that mean this
            # has already been created! Then, skip this step
            if len( total_vars) == 0:
                for i in range( len( input_vars)):
                    temp_list = []
                    # make empty lists for midpoints
                    for j in range( len( midpoints)):
                        temp_list.append( [])
                    total_vars.append( temp_list)

            # save all relevant fields as numpy arrays rather than xarray: saves computing time!
            # also choose whether to sort by rmw or distances in this line of code! keep calling
            # the values "rmwvals" for either case, i'm too lazy to update code names right now lol
            if bintype == 'rmw':
                rmwvals = fl_data['rmw'].values
            elif bintype == 'dist':
                rmwvals = fl_data['center_dist'].values

            datalist = []
            for i in range( len( input_vars)):
                datalist.append( fl_data[ input_vars[i]].values)

            # sort through every rmw axis value for a given eye pass: check if it's within the bin!
            for rmw_i, rmw_val in enumerate( rmwvals ):

                # do nothing for nans
                if np.isnan( rmw_val):
                    nancount += 1
                else:
                    # new method: use array subtraction to find the closest midpoint to the given rmw_i
                    bin_i = np.argmin( np.abs( midpointarray - rmw_val ))

                    # do this for every input variable!
                    for var_i, var_val in enumerate( datalist):
                        # add the singular value for this variable to the correct spot in the total binned array
                        total_vars[ var_i][ bin_i].append( var_val[ rmw_i])

            # transfer the list of lists of lists (lol) to a nice pandas format!
            for col_i, col_val in enumerate( total_vars):
                df_new_bins[ input_vars[ col_i ] ] = col_val

    return df_new_bins

# set up for creating nice mean plots by calculating the mean and confidence limits!
# add them to the dataframe
def plot_setup( df_bin, confidence = .95):
    # do this for all remaining fields
    removekey = ['bins', 'midpoints', 'center_dist', 'rmw']
    templist = set( df_bin.keys() )
    for key in removekey:
        templist.remove( key)
    fieldlist = list( templist )

    for fieldi, fieldval in enumerate( fieldlist):
        var_mean, var_lowc, var_highc = [], [], []
        field = df_bin[ fieldval]

        # do this for every bin for that variable
        for i in range( len( df_bin[ 'midpoints'] )):
            mean, lowc, highc = helpful_stats.t_test_intervals( field[i], confidence=confidence)
            var_mean.append(  mean)
            var_lowc.append(  lowc)
            var_highc.append(  highc)

        # after going through every pass for this variable, add the results to the dataframe!
        df_bin[ fieldval + "_mean"] = var_mean
        df_bin[ fieldval + "_highc"] = var_highc
        df_bin[ fieldval + "_lowc"] = var_lowc

        # print( "Confidence intervals found for " + var_names[ var_i])
    return df_bin
