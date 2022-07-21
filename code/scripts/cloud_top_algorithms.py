import numpy as np
import xarray as xr
from scipy.signal import find_peaks


def cta_top_layer_lowest_value( power, cutoff_power, H_index_matrix, xaxis ):
    # range of power values after these two steps: between -.6 and cutoff dBz
    temp = power.where( power.values > cutoff_power)

    # if there is a nan value, turn it into a zero. Turn non-NaN values into
    # height indices. These height indices will be manipulated in the following code.
    power_threshold = np.where( np.isnan( temp) , 0, H_index_matrix)

    # make an empty array to hold power indices for plotting
    power_index = np.empty( len( xaxis), dtype=int)
    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):
        power_col = power_threshold[ column_index, :]
        # maybe change this, or at least get rid of the fixed number...
        # see selection script for full description
        if power_col[10] == 0:
            power_index[ column_index] = 0

        # remove all zero values
        power_col = power_col[ power_col > 0]

        # check for empty list
        if np.size( power_col) == 0:
            power_index[ column_index] = 0
        else:
            # split current list of indices into arrays of consecutive values! Helps determine the top cloud chunk
            splits = np.split( power_col, np.where(np.diff( power_col) != 1)[0] + 1)
            # look at the last value in the first array of non NaN values: this represents
            # the end of the first clear air chunk beneath the plane! Aka the cloud top
            # add 1 to go from bottom of clear air to top of cloud... actually the
            # + 1 caused an out of bounds error so I got rid of it
            power_index[ column_index] = splits[0] [ -1]
            # power_index[ column_index] = power_col[ -1] # old solution
    return power_index

def cta_lowest_value( power, cutoff_power, H_index_matrix, xaxis):
    # range of power values after these two steps: between -.6 and cutoff dBz
    temp = power.where( power.values > cutoff_power)

    # if there is a nan value, turn it into a zero. Turn non-NaN values into
    # height indices. These height indices will be manipulated in the following code.
    power_threshold = np.where( np.isnan( temp) , 0, H_index_matrix)

    # original cloud top height selection script
    # make an empty array to hold power indices for plotting
    power_index = np.empty( len( xaxis), dtype=int)
    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):
        power_col = power_threshold[ column_index, :]
        power_col = power_col[ power_col > 0]
        # check for empty list
        if np.size( power_col) == 0:
            power_index[ column_index] = 0
        else:
            power_index[ column_index] = power_col[ -1]
    return power_index


def cta_max_value( power, cutoff_power, xaxis):

    power = power.where( power.values > cutoff_power)

    # get rid of nans, they were causing errors for some reason. -200 is far below
    # every max peak
    power = np.where( np.isnan( power) , -200, power)
    power_index = np.empty( len( xaxis), dtype=int)
    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):
        power_col = power[ column_index, :]
        # argmax() returns the index of the max value!
        # print( power_col)
        max_power_ind = np.argmax( power_col)
        power_index[ column_index] = max_power_ind
    return power_index


def cta_find_peaks( power, cutoff_power, xaxis):

    power = power.where( power.values > cutoff_power)
    # get rid of nans, they were causing errors for some reason. -200 is far below
    # every max peak
    power = np.where( np.isnan( power) , -200, power)
    power_index = np.empty( len( xaxis), dtype=int)
    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):
        power_col = power[ column_index, :]

        peaks = find_peaks( power_col, height = -17)
        # pick the last viable power index using [-1]

        # peaks[0] returns an array of indices of values that have peaks!
        # case where no local peaks above -14 dBz exist
        # update this to not just pick the max, maybe use some of the methods above?
        if np.size( peaks[0]) == 0:
            max_power_ind = np.argmax( power_col)
            power_index[ column_index] = max_power_ind
        # case where local peaks exist!
        else:
            # pick the first viable power index using [ 0]
            power_index[ column_index] = peaks[0][0]
    return power_index



def cta_find_peaks_max( power, cutoff_power, xaxis, H_index_matrix):
    power = power.where( power.values > cutoff_power)
    power_threshold = np.where( np.isnan( power) , 0, H_index_matrix)
    # get rid of nans, they were causing errors for some reason. -200 is far below
    # every max peak
    power = np.where( np.isnan( power) , -200, power)
    power_index = np.empty( len( xaxis), dtype=int)

    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):
        power_col = power[ column_index, :]

        # step 1: look for the largest peaks / brightest cloud tops
        peaks = find_peaks( power_col, height = -12)

        # peaks[0] returns an array of indices of values that have peaks!
        # no cloud brighter than -10 dBz:
        if np.size( peaks[0]) == 0:

            # try again for lower power value! Looking for less bright clouds
            # this method is nice because it preferentially chooses the bright clouds
            peaks = find_peaks( power_col, height=-17)

            # if no clouds match this less restrictive threshold, then just choose
            # the last clear sky value from the first continuous dataset = cloud tops
            if np.size( peaks[0] ) == 0:

                # this is just a copy paste of the code used in the cloud_top_lowest_layer() script
                power_col = power_threshold[ column_index, :]
                power_col = power_col[ power_col > 0]
                # check for empty list: no backscattered values at all
                if np.size( power_col) == 0:
                    power_index[ column_index] = 0
                    continue
                # otherwise, pick the lowest backscattered value in the first clear air patch
                else:
                    # split current list of indices into arrays of consecutive values! Helps determine the top cloud chunk
                    splits = np.split( power_col, np.where(np.diff( power_col) != 1)[0] + 1)
                    # look at the last value in the first array of non NaN values: this represents
                    # the end of the first clear air chunk beneath the plane! Aka the cloud top
                    power_index[ column_index] = splits[0] [ -1]
                    continue
            # -10 > power > -16 dBz case
            # pick the first viable power index using [ 0]
            else:
                power_index[ column_index] = peaks[0][0]
                continue
        # -10 < power case
        # pick the first viable power index using [ 0]
        else:
            power_index[ column_index] = peaks[0][0]
            continue

    return power_index





def cta_prominence( power, cutoff_power, xaxis, H_index_matrix):
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
    :return power_index: an array of indices that represent the cloud top heights.
            H[power_index] is called after this function to get the height values.
    """

    # filter the power and power indices from H_index_matrix
    power = power.where( power.values > cutoff_power)
    # this line changes all NaNs to 0s, and all non Nans to indices
    power_ind_matrix = np.where( np.isnan( power) , 0, H_index_matrix)

    # get rid of nans in the power matrix, they were causing errors for some reason.
    # -200 is far below every max peak.
    power = np.where( np.isnan( power) , -200, power)

    # this empty array will be returned at the end of the function call
    power_index = np.empty( len( xaxis), dtype=int)

    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):

        # select the current column from power and power_ind_matrix. Then, get rid of all 0s
        power_column = power[ column_index, :]
        power_ind_column = power_ind_matrix[ column_index, :]
        power_ind_column = power_ind_column[ power_ind_column > 0]

        # check for empty list: no backscattered values at all
        # in this case, the height should have an index of 0 aka the top value
        if np.size( power_ind_column) == 0:
            power_index[ column_index] = 0
            continue

        # otherwise, only look at the first clear air patch
        # split current list of indices into arrays of consecutive values! Helps determine the top cloud chunk
        splits = np.split( power_ind_column, np.where(np.diff( power_ind_column) != 1)[0] + 1)
        # only look at the first clear air patch
        first_consec_inds = splits[0]

        # if the first clear air patch is below heavy attenuation (no signal until 51st height value),
        # just return an index of 0 aka top value
        if first_consec_inds[ 0] > 50:
            power_index[ column_index] = 0
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

            # look at the last value in the first array of non NaN values: this represents
            # the end of the first clear air chunk beneath the plane! Aka the cloud top
            power_index[ column_index] = first_consec_inds[ -1]


        # case where local peaks exist! Use the prominence output from find_peaks
        else:
            # pick the first viable power index using [ 0]
            power_index[ column_index] = peaks[0][0]

    return power_index
