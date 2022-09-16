import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import cloud_top_algorithms as cta
import warnings

# from scipy.signal import find_peaks

def find_cloud_heights( crl_name, cutoff_power, i1, i2, xaxis='time', crl_path = "/Users/etmu9498/research/data/CRL_data/2021", new_heights=False):
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
    warnings.filterwarnings("ignore")

    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)

    if xaxis == 'lon':
        axis = data.Lon[ i1: i2]
    elif xaxis == 'lat':
        axis = data.Lat[ i1: i2]
    elif xaxis == 'time':
        axis = data.time[ i1: i2]
    elif xaxis == 'in-situ-dist':
        axis = data.in_situ_distance[ i1: i2]
    elif xaxis == 'tdr-dist':
        axis = data.tdr_distance[ i1: i2]
    else:
        print( "please select 'lon', 'lat', or 'time' as a valid x axis")

    # new data case: updated height values and corrected power values
    if new_heights:
        H = data.H_new
        power = data.power_new[ i1: i2]
    # old data case: still need to correct power values
    else:
        H = data.H
        power = data.p_ch1[ i1: i2]
        # convert to dBz
        power = 10 * np.log10( power)


    # an index for each of the 594 height values
    H_index = range( len( H))
    # a matrix of height indices, to be used later. 594 height values x 300 repeats
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( axis), axis=0)

    # cut off top 8 values to get rid of NaNs from flight level data

    ''' maybe uncomment these lines! maybe important '''
    # power = power[:, 8:]
    # H_index_matrix = H_index_matrix[:, 8:]

    # power_index = cta.cta_top_layer_lowest_value( power, cutoff_power, H_index_matrix, axis)
    # power_index = cta.cta_max_value( power, cutoff_power, axis)
    # power_index = cta.cta_find_peaks( power, cutoff_power, axis)
    # power_index = cta.cta_find_peaks_max( power, cutoff_power, axis, H_index_matrix)

    '''
    print( len( H))
    print( cutoff_power)
    print( np.shape( power))
    np.set_printoptions(threshold=np.inf)
    print( power.values)
    np.set_printoptions(threshold=1000)
    print( axis.values)
    '''
    if new_heights:
        p3_heights = data.p3_heights[ i1:i2].values / 1000 # convert to km

        # special cases for troublesome grace 8/17 and 8/18 data
        if crl_name[0:15] == 'crl-grace-08-17':
            blank_vals = - p3_heights + np.nanmax( -H) -.4 # this value matches the one in get_p3_heights.py
            grace_case = 1

        elif crl_name[0:15] == 'crl-grace-08-18':
            blank_vals = - p3_heights + np.nanmax( -H) + .65 # this value matches the one in get_p3_heights.py
            grace_case = 2

        # normal case
        else:
            grace_case = False
            # blank_vals should be close to - 3.25 + 4.0 = .75
            blank_vals = - p3_heights + np.nanmax( -H)


        power_index = cta.cta_prominence_in_situ( power, cutoff_power, axis, H_index_matrix, p3_heights, -H, grace_case)

        # the heights returned by the power indices need to be corrected! Before height
        # corrections were made, an index of 0 would automatically correspond to a height of
        # 3.557 km, which worked! But now, theres a ton of blank space from 4 km down to flight
        # level. So, we can no longer assume that an index of 0 should be at the top, aka 4 km.

        # This code corrects for this height difference by finding the width of the blank
        # space and subtracting it from the returned height to shift the data down.


        cloud_heights = -H[ power_index] - blank_vals
        cloud_heights = np.where( cloud_heights < 0, 0, cloud_heights)
    else:
        power_index = ta.cta_prominence( power, cutoff_power, axis, H_index_matrix)
        cloud_heights = -H[ power_index]

    '''
    print( power_index)
    print( - H[ power_index].values)

    print( "\n" + str( len( H)))
    print( H.values)
    '''

    warnings.filterwarnings("default")

    return cloud_heights, axis





# the same as the function above, but it finds multiple cloud layers, not just the highest!
def find_multi_cloud_heights( crl_name, cutoff_power, i1, i2, xaxis='in-situ-dist', crl_path = "/Users/etmu9498/research/data/crl-new-matrices"):
    warnings.filterwarnings("ignore")
    # load data
    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)
    # load x axis data
    if xaxis == 'in-situ-dist':
        axis = data.in_situ_distance[ i1: i2]
    else:
        print( "update find_multi_cloud_heights() if statement!")

    # new data case: updated height values and corrected power values
    H = data.H_new
    power = data.power_new[ i1: i2]

    # an index for each of the height values
    H_index = range( len( H))
    # a matrix of height indices, to be used in cloud top algorithm helper fn
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( axis), axis=0)
    # get p3 heights
    p3_heights = data.p3_heights[ i1:i2].values / 1000 # convert to km

    # account for 'blank values', aka the nans below the p-3 height, to properly scale
    # the cloud top heights.
    # special case for troublesome grace 8/18 data: this dataset was shifted by .65 km,
    # so we also need to shift height values by the same amount
    if crl_name[0:15] == 'crl-grace-08-18':
        blank_vals = - p3_heights + np.nanmax( -H) + .65 # this value matches the one in get_p3_heights.py
        grace_case = True
    # normal case
    else:
        grace_case = False
        # blank_vals should be close to - 3.25 + 4.0 = .75
        blank_vals = - p3_heights + np.nanmax( -H)

    # call the helper function to get cloud height locations
    power_index = cta.cta_prominence_many_cloud_layers( power, cutoff_power, axis, H_index_matrix, p3_heights, -H, grace_case)

    # initialize variables for returns later!
    cloud_heights = [] # height of clouds at each x axis value
    new_xaxis = [] # x axis values, with repeating values for multiple cloud height cases
    cloud_counts = np.zeros( 100).tolist() # number of cases with 0, 1, 2, etc cloud layers

    # Find cloud heights and x axis values for each index... the for loop is needed
    # to account for multiple values at each layer!
    for i in range( len( power_index)):
        cloud_heights += [ -H[ power_index[i] ].values - blank_vals[ i] ]
        new_xaxis += np.full( shape=len(power_index[i]), fill_value= axis[i] ).tolist()

        # also, count the numebr of cloud layers!
        cloud_number_i = len( power_index[ i])
        cloud_counts[ cloud_number_i] += 1

    # trim the huge number of trailing zeros off the array!
    cloud_counts = np.trim_zeros( np.array( cloud_counts), trim='b').tolist()

    # try flattening all of the lists for easier plotting
    # this turns things from a list of lists into one flattened list!
    flat_power_index = [item for sublist in power_index for item in sublist]
    flat_cloud_heights = [item for sublist in cloud_heights for item in sublist]

    warnings.filterwarnings("default")
    return flat_cloud_heights, new_xaxis, cloud_counts



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
        leg = plt.legend( loc='lower left')
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
        leg = plt.legend( loc='lower left')
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
