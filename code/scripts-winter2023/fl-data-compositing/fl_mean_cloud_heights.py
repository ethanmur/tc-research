import numpy as np
import xarray as xr
import os
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import scipy
import sys
import warnings
import pandas as pd
from scipy.signal import find_peaks

sys.path.append("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns
import helpful_stats
os.chdir("/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import find_cloud_tops

# like other flight level datasets, calculate the mean cloud height vs rmw!
# use the crl data instead of flight level data
def calc_cloud_heights_rmws( tc='all', binwidth=.1, maxbin=10, confidence=.95):
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
    cloud_heights_df = []
    p3_heights_df = []
    # initial case: make a bunch of empty lists in cloud_heights_df and p3_heights_df!
    # will be filled below in for loops.
    # make empty lists for midpoints
    for j in range( len( midpoints)):
        cloud_heights_df.append([])
        p3_heights_df.append([])

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            
            # get data
            crl_path = crl_data_root + yearval
            os.chdir( crl_path)
            crl_data = xr.open_dataset( fileval)

            # save all relevant fields as numpy arrays rather than xarray: saves computing time!
            rmwvals = crl_data['rmw'].values
            p3_heightvals = crl_data['p3_height'].values

            # find cloud heights for this case!
            H = crl_data.height
            power = crl_data.P_ch1
            axis = crl_data.time
            p3_height = crl_data.p3_height
            if yearval == '2021':
                min = -30
            elif yearval == '2022':
                min = -40
            cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = min)

            print("Cloud heights found")

            # sort through every rmw axis value for a given eye pass: check if it's within the bin!
            for rmw_i, rmw_val in enumerate( rmwvals ):

                # do nothing for nans
                if np.isnan( rmw_val):
                    nancount += 1
                else:
                    # new method: use array subtraction to find the closest midpoint to the given rmw_i
                    bin_i = np.argmin( np.abs( midpointarray - rmw_val ))

                    # add the singular value for this variable to the correct spot in the total binned array
                    cloud_heights_df[ bin_i].append( cloudheights[ rmw_i])

                    # add the p3 height at this point to a similar array
                    p3_heights_df[ bin_i].append( p3_heightvals[ rmw_i])

    df_new_bins['cloud_heights'] = cloud_heights_df
    df_new_bins['p3_heights'] = p3_heights_df

    print( "Cloud and P-3 Heights Added to Dataframe")

    # calculate confidence intervals for the new cloud height variables
    var_mean, var_lowc, var_highc = [], [], []
    p3_mean = []
    field = df_new_bins[ 'cloud_heights']

    # do this for every bin for that variable
    for i in range( len( df_new_bins[ 'midpoints'] )):
        mean, lowc, highc = helpful_stats.t_test_intervals( field[i], confidence=confidence)
        var_mean.append(  mean)
        var_lowc.append(  lowc)
        var_highc.append(  highc)
        p3_mean.append( np.mean(df_new_bins['p3_heights'][i]))

    # after going through every pass for this variable, add the results to the dataframe!
    df_new_bins["cloud_heights_mean"] = var_mean
    df_new_bins["cloud_heights_highc"] = var_highc
    df_new_bins["cloud_heights_lowc"] = var_lowc

    df_new_bins['p3_heights_mean'] = p3_mean

    print( "Means and Confidence Intervals Found")

    return df_new_bins






# like the code above, but calculate these cloud height distributions for individual passes!
# this might help highlight
# use given eyewall limits for rmw = 1 locations?
def calc_cloud_heights_rmws_single_passes( tc='all', binwidth=.1, maxbin=10, confidence=.95):
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
    cloud_heights_df = []
    # initial case: make a bunch of empty lists in cloud_heights_df
    # make empty lists for midpoints
    for j in range( len( midpoints)):
        cloud_heights_df.append([])

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            
            print(fileval)

            # get data
            crl_path = crl_data_root + yearval
            os.chdir( crl_path)
            crl_data = xr.open_dataset( fileval)

            # save all relevant fields as numpy arrays rather than xarray: saves computing time!
            rmwvals = crl_data['rmw'].values
            # find cloud heights for this case!
            H = crl_data.height
            power = crl_data.P_ch1
            axis = crl_data.time
            p3_height = crl_data.p3_height
            if yearval == '2021':
                min = -30
            elif yearval == '2022':
                min = -40
            cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = min)

            # sort through every rmw axis value for a given eye pass: check if it's within the bin!
            for rmw_i, rmw_val in enumerate( rmwvals ):

                # do nothing for nans
                if np.isnan( rmw_val):
                    nancount += 1
                else:
                    # new method: use array subtraction to find the closest midpoint to the given rmw_i
                    bin_i = np.argmin( np.abs( midpointarray - rmw_val ))

                    # add the singular value for this variable to the correct spot in the total binned array
                    cloud_heights_df[ bin_i].append( cloudheights[ rmw_i])

    df_new_bins['cloud_heights'] = cloud_heights_df

    # calculate confidence intervals for the new cloud height variables
    var_mean, var_lowc, var_highc = [], [], []
    field = df_new_bins[ 'cloud_heights']

    # do this for every bin for that variable
    for i in range( len( df_new_bins[ 'midpoints'] )):
        mean, lowc, highc = helpful_stats.t_test_intervals( field[i], confidence=confidence)
        var_mean.append(  mean)
        var_lowc.append(  lowc)
        var_highc.append(  highc)

    # after going through every pass for this variable, add the results to the dataframe!
    df_new_bins["cloud_heights_mean"] = var_mean
    df_new_bins["cloud_heights_highc"] = var_highc
    df_new_bins["cloud_heights_lowc"] = var_lowc


    return df_new_bins