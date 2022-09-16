"""
Created on Mon Aug  15 13:58 2022
create and return an array of p-3 heights for every crl x axis value
@author: ethanmur
"""
import os
import xarray as xr
import numpy as np
os.chdir( "/Users/etmu9498/research/code/scripts/")
import tc_metadata
os.chdir( "/Users/etmu9498/research/code/scripts/in-situ-scripts")
import load_in_situ_data

def geth( tcname, dataset):
    metadata = tc_metadata.all_data( tc= tcname)
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    in_situ_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

    # load data
    in_situ_path = "/Users/etmu9498/research/data/in-situ-new"
    crl_path = "/Users/etmu9498/research/data/crl-new"
    os.chdir( in_situ_path)
    in_situ_data = xr.open_dataset( in_situ_name)
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # in situ time in hours
    is_time_hours = in_situ_data.float_time
    # time in seconds (UTC)
    is_time_secs = in_situ_data.float_time.values * 60 * 60

    # define the height of the P-3 above ocean surface
    # height = [ float( line) for line in in_situ_data["HT.d"].values ]
    # not as elegant as above. workaround for dealing with empty strings
    height = in_situ_data["HT.d"].values
    htemp = np.zeros( len( height))
    for line_ind in range( len( height)):
        if height[ line_ind] == '':
            htemp[ line_ind] = np.nan
        else:
            htemp[ line_ind] = float( height[ line_ind])
    height = htemp.tolist()

    # make an empty list to store p-3 height values
    # this will be the length of the crl time, lat, etc variables
    height_list = []
    # do this for every time step in the crl data
    for i in range( len( crl_data.time)):

        # find the index of the closest in situ time value to the current crl time
        index = (np.abs( is_time_hours.values - crl_data.time[ i].values)).argmin()

        # find the height value at that corresponding time
        # add it to the list
        height_i = height[ index]
        height_list += [ height_i]

    return height_list

# this function interprets the data in the following way:
# it creates a new height matrix from 0 to the p_3 height and 'squishes'
# the data into this smaller array
# resolution is in meters here
def interp_data(matrix, p3_heights, resolution=6):
    # defining a new base height array to standardize the new matrix
    base_h = np.arange( -4.0, 0.0, resolution / 1000)
    # make an empty array for initial storage
    new_matrix=np.empty([ np.size( matrix, 0), len( base_h)])

    # do this for every x axis value
    for i in range( np.size( matrix, 0) ):
        # get the current matrix column
        columni = matrix[ i, :]

        # make a new height profile for this column, from 0 to p-3 height
        trim_h = np.linspace( - p3_heights[ i] / 1000, 0.0, np.size( matrix, 1))

        # interpolate over that row and save results
        new_columni = np.interp( base_h, trim_h, columni)
        new_matrix[ i, :] = new_columni

    return base_h, new_matrix

# this function interprets the data in the following way:
# It finds the difference between the original crl height ( ~3.45 km) and the
# flight height ( ~ 3.0 km). This difference should account for the increased bottom
# height in the rainbands.
# instead of 'squishing' the data, this function shifts the data down by the
# difference that was just found
# resolution is in meters here

# temp max isn't important; just another test to fix bottom crl height limit
def interp_data2(matrix, original_crl_heights, p3_heights, grace_case=False, temp_max=None, resolution=6):
    # defining a new base height array to standardize the new matrix
    base_h = np.arange( -4.0, 0.0, resolution / 1000)
    # make an empty array for initial storage
    new_matrix=np.empty([ np.size( matrix, 0), len( base_h)])

    # np.set_printoptions(threshold=np.inf)
    # print( matrix[ 0:10, :].values)
    # np.set_printoptions(threshold=1000)


    # find the max height value with actual data for the crl data
    # this should be the same for every tc case
    # orig_maxh = np.nanmax( - original_crl_heights)

    # the following code is actually all unnecesary...

    # redefine variable
    origh = original_crl_heights
    orig_maxh = np.nanmax( - origh)

    '''
    # the top few values are filled with 999. for some reason. find the maxh
    # without 999. values present. Actually, they're filled with nans: sort those

    # look at x axis column 10 or 200 to get rid of nans effectively; maybe change this number
    # if things aren't working!
    # case where issues at column 10 arise
    if len( np.where( ~ np.isnan (matrix[ 10, :] ) )[0] ) == 0:
        orig_maxh = np.nanmax( - origh[  np.where( ~ np.isnan (matrix[ 200, :] ) ) ])
    # base case: this is most likely
    else:
        orig_maxh = np.nanmax( - origh[  np.where( ~ np.isnan (matrix[ 10, :] ) ) ])
    '''

    # do this for every x axis value
    for i in range( np.size( matrix, 0) ):
        # get the current matrix column
        columni = matrix[ i, :]

        # calculate a couple more useful heights: height from p3 and difference
        # from original estimate
        p3h = p3_heights[ i] / 1000
        shift = 0.0
        # when the top Nans in the T array ARE included (the right way, I think),
        # a 50 m vertical shift is needed to reproduce the results below.
        # a 75 m shift seems to keep some of the bright reflectivity scattering
        # I'm still gonna assume 0 m vertical shift for now... to not mess with results?

        # when the top 25 or so nans AREN'T included:
        # it seems like a shift of .05 km includes the strong reflectivity
        # a shift of .010 (only 10 m!) gets rid of most the reflectivity
        # and no shift might get rid of 10 m of valid data, but is the best
        # method I can think of for trimming data

        # need to shift data down by an extra 400 m
        if grace_case == 17:
            hdiff = orig_maxh - p3h - shift - .4 # old shift was .95
            trim_h = np.linspace( - p3h - .4 , 0.0 + hdiff, np.size( matrix, 1))

        # need to shift data up by an extra 650 m
        elif grace_case == 18:
            hdiff = orig_maxh - p3h - shift + .65 # old shift was .95
            trim_h = np.linspace( - p3h + .65 , 0.0 + hdiff, np.size( matrix, 1))

        else:
            hdiff = orig_maxh - p3h - shift # the - shift is a lil extra padding to see what's going on
            # make a new height profile for this column, from 0 to p-3 height
            # this is different from the step above!
            # ex: - p3h ~ -3.15, 0 - hdiff ~ + .4, should get rid of lower vals!
            trim_h = np.linspace( - p3h , 0.0 + hdiff, np.size( matrix, 1))


        '''
        print( 'step: ' + str( i))
        print( 'original p3 height: ' + str( orig_maxh))
        print( 'new p3 height: ' + str( p3h))
        print( hdiff)
        '''


        # interpolate over that row and save results
        new_columni = np.interp( base_h, trim_h, columni)
        new_matrix[ i, :] = new_columni

    return base_h, new_matrix
