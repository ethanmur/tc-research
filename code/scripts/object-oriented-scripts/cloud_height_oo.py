import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
from scipy.signal import find_peaks



def cloud_top_alg_top_layer_lowest_value( power, cutoff_power, H_index_matrix, xaxis ):
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

def cloud_top_alg_lowest_value( power, cutoff_power, H_index_matrix, xaxis):
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


def cloud_top_alg_max_value( power, cutoff_power, xaxis):

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


def cloud_top_alg_find_peaks( power, cutoff_power, xaxis):

    power = power.where( power.values > cutoff_power)
    # get rid of nans, they were causing errors for some reason. -200 is far below
    # every max peak
    power = np.where( np.isnan( power) , -200, power)
    power_index = np.empty( len( xaxis), dtype=int)
    # cycle through each column until the last power value within cutoff is found: this is the cloud top height!
    for column_index in range( len( xaxis)):
        power_col = power[ column_index, :]

        peaks = find_peaks( power_col, height = -14)
        # pick the last viable power index using [-1]

        # peaks[0] returns an array of indices of values that have peaks!
        # case where no local peaks above -14 dBz exist
        # update this to not just pick the max, maybe use some of the methods above?
        if np.size( peaks[0]) == 0:
            max_power_ind = np.argmax( power_col)
            power_index[ column_index] = max_power_ind
        # case where local peaks exist!
        else:
            power_index[ column_index] = peaks[0][-1]
    return power_index


def find_cloud_heights( crl_name, cutoff_power, i1, i2, xaxis='time'):
    """
    This function finds and returns cloud heights
    :param crl_name: The name of the crl data file that will be analyzed. The crl_path
            is assumed to remain fixed for all datasets
    :param cutoff_power: The minimum threshold power that the backscattered signal
            must meet to be included in the plot. Anything smaller than this number will be removed.
    :return power_time: A list of times from the eyewall pass that will be used for plotting.
    :return power_height: A list of heights representing the lowest layer of clear air detected
            by the lidar
    """
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)
    power = data.P_ch1[ i1: i2]

    if xaxis == 'lon':
        axis = data.Lon[ i1: i2]
    elif xaxis == 'lat':
        axis = data.Lat[ i1: i2]
    elif xaxis == 'time':
        axis = data.time[ i1: i2]
    else:
        print( "please select 'lon', 'lat', or 'time' as a valid x axis")

    H = data.H
    # an index for each of the 594 height values
    H_index = range( len( H))
    # a matrix of height indices, to be used later. 594 height values x 300 repeats
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( axis), axis=0)

    # cut off top 8 values to get rid of NaNs from flight level data
    power = power[:, 8:]
    H_index_matrix = H_index_matrix[:, 8:]

    # convert to dBz
    power = 10 * np.log10( power)

    # power_index = cloud_top_alg_top_layer_lowest_value( power, cutoff_power, H_index_matrix, axis)
    power_index = cloud_top_alg_max_value( power, cutoff_power, axis)
    # power_index = cloud_top_alg_find_peaks( power, cutoff_power, axis)
    return - H[power_index], axis



def cloud_height_crl_comparison( crl_name, cutoff_power, i1, i2, xaxis='time', xlims=False, same_plot=False):
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"

    H, xaxis_value = find_cloud_heights(crl_name, cutoff_power, i1, i2, xaxis = xaxis)
    plt.clf()

    # plot cloud top height line right on top of crl data
    line_color = 'g'
    if same_plot:
        # plot results
        plt.figure( figsize=(18, 5))
        make_plots.plot_power_ch1( crl_path, crl_name, i1, i2, xaxis, cutoff= cutoff_power)
        plt.scatter( xaxis_value, H, c= line_color, s=8, marker='s') # s
        plt.plot( xaxis_value, H, c=  line_color, linewidth=2, label= 'Cloud Top Height')
        if xlims:
            plt.xlim( [ xlims[0], xlims[1] ])
        # make lines in legend thicker
        leg = plt.legend()
        for line in leg.get_lines():
            line.set_linewidth(4.0)
    # plot cloud top height data and crl data separately
    else:
        # plot results
        plt.figure( figsize=(14.9, 4))
        plt.scatter( xaxis_value, H, c= line_color, s=12, marker='s')
        plt.plot( xaxis_value, H, c= line_color, linewidth=.5, label= 'Cloud Top Height')
        if xlims:
            plt.xlim( [ xlims[0], xlims[1] ])
        plt.ylabel( "height (km)")
        plt.xlabel( 'time (hours)')
        plt.grid('on')
        # make lines in legend thicker
        leg = plt.legend()
        for line in leg.get_lines():
            line.set_linewidth(4.0)

        plt.figure( figsize=(29, 4))
        make_plots.plot_power_ch1( crl_path, crl_name, i1, i2, xaxis, cutoff= cutoff_power)
        if xlims:
            plt.xlim( [ xlims[0], xlims[1] ])




def cloud_height_sensitivity( crl_name, cutoff_list, xlims=False, colors=False, labels=False):
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)
    data

    lon = data.Lon
    lat = data.Lat
    power = data.P_ch1
    time = data.time
    H = data.H
    # an index for each of the 594 height values
    H_index = range( len( H))
    # a matrix of height indices, to be used later. 594 height values x 300 repeats
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( time), axis=0)


    # cut off top 8 values to get rid of error from flight level data
    power = power[:, 8:]
    H_index_matrix = H_index_matrix[:, 8:]

    # convert to dBz
    power = 10 * np.log10( power)

    plt.figure( figsize=(14.9, 4))
    for i in range( len( cutoff_list)):
        power_cutoff = cutoff_list[ i]
        # power_index = cloud_top_alg_top_layer_lowest_value( power, power_cutoff, H_index_matrix, time)
        power_index = cloud_top_alg_max_value( power, power_cutoff, time)
        # plot results
        if not colors and not labels:
            plt.scatter( time, - H[ power_index], s=8, marker='s')
            plt.plot( time, - H[power_index], linewidth=.5, label= cutoff_list[ i])
        elif not colors:
            plt.scatter( time, - H[ power_index], s=8, marker='s')
            plt.plot( time, - H[power_index], linewidth=.5, label= labels[ i])
        elif not labels:
            plt.scatter( time, - H[ power_index], c = colors[i], s=8, marker='s')
            plt.plot( time, - H[power_index], c = colors[i], linewidth=.5, label= cutoff_list[ i])
        else:
            plt.scatter( time, - H[ power_index], c = colors[i], s=8, marker='s')
            plt.plot( time, - H[power_index], c = colors[i], linewidth=.5, label= labels[ i])

    plt.grid('on')
    plt.ylabel( 'Height (km)')
    plt.xlabel( 'Time (hours)')
    if xlims:
        plt.xlim( [ xlims[0], xlims[1] ])

    # make lines in legend thicker
    leg = plt.legend()
    for line in leg.get_lines():
        line.set_linewidth(4.0)






def algorithm_comparison( crl_name, cutoff_power):
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)
    data

    lon = data.Lon
    lat = data.Lat
    power = data.P_ch1
    time = data.time
    H = data.H
    # an index for each of the 594 height values
    H_index = range( len( H))
    # a matrix of height indices, to be used later. 594 height values x 300 repeats
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( time), axis=0)

    # cut off top 8 values to get rid of NaNs from flight level data
    power = power[:, 8:]
    H_index_matrix = H_index_matrix[:, 8:]
    # convert to dBz
    power = 10 * np.log10( power)

    # power_index = cloud_top_alg_top_layer_lowest_value( power, cutoff_power, H_index_matrix, time)
    power_index_max = cloud_top_alg_max_value( power, cutoff_power, time)
    power_index_top_layer = cloud_top_alg_top_layer_lowest_value( power, cutoff_power, H_index_matrix, time )
    power_index_lowest_layer = cloud_top_alg_lowest_value( power, cutoff_power, H_index_matrix, time)
    power_index_find_peaks = cloud_top_alg_find_peaks( power, cutoff_power, time)

    return time, - H[power_index_max], - H[power_index_top_layer], - H[power_index_lowest_layer], - H[power_index_find_peaks]
