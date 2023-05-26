import os
import numpy as np
import xarray as xr
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts/rainfall")
import rainfall_algorithm
os.chdir(  "/Users/etmu9498/research/code/scripts")
import cloud_height
import cloud_top_algorithms as cta



# this function finds regions of rainfall for plotting and statisical purposes
def find_rain( crl_name, i1, i2, cutoff_power=-30, xaxis='in-situ-dist', crl_path = "/Users/etmu9498/research/data/crl-new-matrices"):
    warnings.filterwarnings("ignore")
    # load data
    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)
    # load x axis data
    if xaxis == 'in-situ-dist':
        axis = data.in_situ_distance[ i1:i2]
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
    p3_heights = data.p3_heights[i1:i2].values / 1000 # convert to km

    # account for 'blank values', aka the nans below the p-3 height, to properly scale
    # the cloud top heights.
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

    # call the helper function to get cloud height locations
    cloud_height_inds = cta.cta_prominence_many_cloud_layers( power, cutoff_power, axis, H_index_matrix, p3_heights, -H, grace_case)

    power_index = rainfall_algorithm.find_rain( power, cutoff_power, axis, H_index_matrix, p3_heights, -H, grace_case, cloud_heights=cloud_height_inds)

    # initialize variables for returns later!
    rainfall_points = [] # height of clouds at each x axis value
    new_xaxis = [] # x axis values, with repeating values for multiple cloud height cases


    ''' fix this output code! rainfall_algorithm looks good, things just aren't being properly flattened '''

    # print( axis)
    print( len( power_index))
    print( len( axis))

    # Find cloud heights and x axis values for each index... the for loop is needed
    # to account for multiple values at each layer!
    for i in range( len( power_index)):
        rainfall_points += [ ( -H[ power_index[i] ].values - blank_vals[ i] ).tolist() ]
        new_xaxis += [ np.full( shape=len(power_index[i]), fill_value= axis[i] ).tolist() ]


        '''
        print( len( rainfall_points))
        print( len( new_xaxis))
        print( rainfall_points)
        print( new_xaxis)
        '''

    # try flattening all of the lists for easier plotting
    # this turns things from a list of lists into one flattened list!
    flat_power_index = [item for sublist in power_index for item in sublist]
    flat_rainfall = [item for sublist in rainfall_points for item in sublist]

    flat_new_xaxis = [item for sublist in new_xaxis for item in sublist]

    warnings.filterwarnings("default")
    return flat_rainfall, flat_new_xaxis


# this function finds regions of rainfall for plotting and statisical purposes
def rain_peaks_wrapper( crl_name, i1, i2, cutoff_power=-30, xaxis='in-situ-dist', crl_path = "/Users/etmu9498/research/data/crl-new-matrices"):
    warnings.filterwarnings("ignore")
    # load data
    os.chdir( crl_path)
    data = xr.open_dataset( crl_name)
    # load x axis data
    if xaxis == 'in-situ-dist':
        axis = data.in_situ_distance[ i1:i2]
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
    p3_heights = data.p3_heights[i1:i2].values / 1000 # convert to km

    # account for 'blank values', aka the nans below the p-3 height, to properly scale
    # the cloud top heights.
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

    # call the helper function to get cloud height locations
    cloud_height_inds = cta.cta_prominence_many_cloud_layers( power, cutoff_power, axis, H_index_matrix, p3_heights, -H, grace_case)

    rainfall_index = rainfall_algorithm.find_rain_peaks( power, cutoff_power, axis, H_index_matrix, p3_heights, -H, grace_case, cloud_heights=cloud_height_inds)

    # initialize variables for returns later!
    rainfall_points = [] # height of clouds at each x axis value
    new_xaxis = [] # x axis values, with repeating values for multiple cloud height cases


    # Find cloud heights and x axis values for each index... the for loop is needed
    # to account for multiple values at each layer!
    for i in range( len( rainfall_index)):
        rainfall_points += [ ( -H[ rainfall_index[i] ].values - blank_vals[ i] ).tolist() ]
        new_xaxis += [ np.full( shape=len( rainfall_index[i]), fill_value= axis[i] ).tolist() ]


    # try flattening all of the lists for easier plotting
    # this turns things from a list of lists into one flattened list!
    flat_power_index = [item for sublist in rainfall_index for item in sublist]
    flat_rainfall = [item for sublist in rainfall_points for item in sublist]

    flat_new_xaxis = [item for sublist in new_xaxis for item in sublist]

    warnings.filterwarnings("default")
    return flat_rainfall, flat_new_xaxis
