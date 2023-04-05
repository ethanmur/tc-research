# created 2/28/23
# code taken from "scripts/cloud_height.py" and

import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.signal import find_peaks

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import cloud_top_algorithms as cta
import warnings


def find_cloud_heights( H, power, axis, p3_height, cutoff_power ):
    """
    This function finds and returns cloud heights
    :param cutoff_power:
    :return power_time:
    :return power_height:
    """
    warnings.filterwarnings("ignore")

    # loading old data
    # H = crl_data.height
    # power = crl_data.P_ch1
    # axis = crl_data.time.values

    # an index for each of the height values
    H_index = range( len( H))
    # a matrix of height indices, to be used later. number of height values x number of repeats
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( axis), axis=0)

    power_index = cloud_top_prominence( power, cutoff_power, axis.values, H_index_matrix, p3_height.values, H)

    matrix_top_to_p3 = -1 * p3_height.values + np.nanmax( H)
    cloud_heights = H[ power_index].values - matrix_top_to_p3
    # make sure there are no negative values
    cloud_heights = np.where( cloud_heights < 0, 0, cloud_heights)

    # cloud_heights = H[ power_index].values

    warnings.filterwarnings("default")
    return cloud_heights, axis




def cloud_top_prominence( power, cutoff_power, xaxis, H_index_matrix, p3_height, Harray):
    """
    This algorithm uses the find_peak() function, with a focus on the prominence
    parameter, to determine
    :param power: A 2D array of backscattered CRL power values
    :param cutoff_power: The minimum threshold power that the backscattered signal
            must meet to be included in the plot. Anything smaller than this number will be removed.
    :param xaxis: A generic name for the x axis values. They can be 'time', 'lat',
            or 'lon', depending on the selection made in the parent function.
    :param H_index_matrix: A 2D array the same size as power. Instead of holding power
            values, though, they hold indices for power values. These indices are
            used to determine gaps in power measurements. The 2D array was made in
            the parent function.
    :param p3_height: the height of the p-3 aircraft. this height is used as the default cloud top height
    :param Harray:
    :return power_index: an array of indices that represent the cloud top heights.
            H[power_index] is called after this function to get the height values.
    """

    # change this to print out more or fewer values!
    # make > 7000 ish to print nothing!
    printind = 10000

    # filter the power and power indices from H_index_matrix
    power = power.where( power.values > cutoff_power)
    # this line changes all NaNs to 0s, and all non Nans to indices
    power_ind_matrix = np.where( np.isnan( power) , 0, H_index_matrix)
    # get rid of nans in the power matrix, they were causing errors for some reason.
    # -200 is far below every max peak.
    power = np.where( np.isnan( power) , -200, power)


    # this empty array will be returned at the end of the function call with the height indices
    power_index = np.empty( len( xaxis), dtype=int)
    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):

        # a very useful helper quantity
        # return the height index closest to this vert profile's p3 height
        Hind = np.argmin( np.abs(Harray.values - p3_height[ column_index]))

        # select the current column from power and power_ind_matrix. Then, get rid of all 0s
        power_column = power[ column_index, :]
        power_ind_column = power_ind_matrix[ column_index, :]
        power_ind_column = power_ind_column[ power_ind_column > 0]

        # check for empty list: no backscattered values at all
        # in this case, the height should have an index of 0 aka the top value
        if np.size( power_ind_column) == 0:
            # just make the cloud height the height of the plane
            power_index[ column_index] = 0 # Hind
            if column_index % printind == 0 and column_index > 0:
                #print( p3_height[ column_index])
                #print( Harray[ Hind])
                print( "Nan Case!\n")
            continue

        # look at the first clear air patch
        # split current list of indices into arrays of consecutive values! Helps determine the top cloud chunk
        splits = np.split( power_ind_column, np.where(np.diff( power_ind_column) != 1)[0] + 1)
        # only look at the first clear air patch
        first_consec_inds = splits[0]
        if column_index % printind == 0 and column_index > 0:
            print( first_consec_inds)


        # similar code as above, but more nuanced! look 50 values below the p-3 height.
        # if it's still a blank section, return the p-3 flight height level!
        if first_consec_inds[0] > Hind + 50:
            # just call the cloud height the height of the plane
            power_index[ column_index] = 0 # Hind
            if column_index % printind == 0 and column_index > 0:
                print( "height correction case\n")
            continue

        # trim the power column to the size of the first clear air patch
        power_column = power_column[ first_consec_inds[0] : first_consec_inds[-1] ]

        # try to find peaks in this dataset!

        # prominence tells the function to look for peaks based on valleys around it.
        # this is good at cutting through to actual peaks! see saved stack overflow post.
        # wlen tells find_peaks to only choose values within a given range for prominence calculations.
        # this gets rid of rain backscatter, as those signals have wide but gradual prominence profiles.
        # the minimum height value provided makes sure that very weak signals aren't included
        # (basically just another clear air filter)

        # width makes sure that there are nearby, tall points around the peak.
        # this gets rid of aerosol backscatter, as these returns are extremely sharp / narrow.
        # 7/21/22 update: the short but sharp peaks are ok, so the width parameter has been removed
        peaks = find_peaks( power_column, prominence= 5, wlen=100, height=-25) # , width = 5)


        # if no peaks exist, pick the end of the first clear cloud layer
        if np.size( peaks[0]) == 0:

            if column_index % printind == 0 and column_index > 0:
                print( 'no peak at i = ' + str( column_index) + ' or xaxis = ' + str( xaxis[ column_index]) + "\n")

            # look at the last value in the first array of non NaN values: this represents
            # the end of the first clear air chunk beneath the plane! Aka the cloud top
            # this value actually needs to get shifted up by the distance of the blank
            # space! This line of code was causing some problems with the new in situ data
            ind_first_power_val = ( np.abs( np.subtract( Harray, p3_height[ column_index]))).argmin()

            # look at the last value in the first array of non NaN values: this represents
            # the end of the first clear air chunk beneath the plane! Aka the cloud top
            power_index[ column_index] = first_consec_inds[ -1] - ind_first_power_val

        # case where local peaks exist! Use the prominence output from find_peaks
        else:
            # pick the first viable power index using [ 0]
            peak = peaks[0][0]
            # NEW
            # shift the peak from 0 down to the height of the p-3! or else points will be too high
            #if first_consec_inds[ -1] + Hind <= len( Harray) - 1:
            power_index[ column_index] = peak # + Hind
            #else:
            #    power_index[ column_index] = len( Harray) - 1


            if column_index % printind == 0 and column_index > 0:
                print( "peak found!\n")

    return power_index
