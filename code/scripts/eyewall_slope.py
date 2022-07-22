import os
import xarray as xr
import numpy as np
from scipy.signal import find_peaks

os.chdir( "/Users/etmu9498/research/code/scripts/")
import make_plots
import plot_in_situ

import warnings



# algorithm for finding points along the eyewall!
def eyewall_slope_first_val( tdr_path, inbound_name, outbound_name, eyewall_cutoff=True, xaxis='dist'):

    # using a cutoff of 20 dBz because it seems like other researchers have used this value for eyewall slope measurements!
    cutoff = 20.0

    # load data
    os.chdir( tdr_path)
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # set up useful variables
    # choose a valid x axis
    if xaxis == 'dist':
        x_in = inbound_data.radius
        x_out = - outbound_data.radius
    elif xaxis == 'lon':

        # code taken from make_plots.py file, plot_tdr() function
        x_in =  inbound_data.longitude [ ~np.isnan( inbound_data.longitude)]
        x_out = outbound_data.longitude [ ~np.isnan( outbound_data.longitude)]

    elif xaxis == 'lat':
        # code taken from make_plots.py file, plot_tdr() function
        x_in =  inbound_data.latitude [ ~np.isnan( inbound_data.latitude)]
        x_out = outbound_data.latitude [ ~np.isnan( outbound_data.latitude)]
    else:
        print( 'choose \'lat\', \'lon\', or \'dist\' for an x axis! ')
        return

    H = inbound_data.height
    H_index = range( len( H))
    # a matrix of height indices, to be used later. number of height values times number of x values
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( x_in), axis=0)


    # work on inbound data first
    in_refl = inbound_data.REFLECTIVITY.isel(time=0).isel(heading=0) # .transpose()

    # only for lat lon cases ( very annoying :/)
    if xaxis == 'lat' or xaxis == 'lon':
        # get rid of nan values from end of columns (like in plot_tdr() script in make_plots.py)
        in_refl = in_refl[ range( len( x_in) ), :]

    # turn all nans to zeros
    in_refl = np.where( np.isnan( in_refl.values), 0, in_refl)

    # return *indices* where our reflectivity value is greater than 20. other values are 0
    in_refl_ind = np.where( in_refl <= cutoff, 0, H_index_matrix )

    # make an empty array to hold eyewall location values
    in_peaks = np.empty( len( x_in), dtype=int)
    # cycle through each column until the first reflectivity value within cutoff is found: this is the start of the eyewall!
    for column_index in range( len( x_in)):
        refl_col = in_refl_ind[ column_index, :]

        # remove zeros: leave only indices with valid reflectivity returns
        refl_col = refl_col[ refl_col > 0]

        # case 1: no dBz values are high enough to meet threshold
        if np.size( refl_col) == 0:
            in_peaks[ column_index] = 0 # float('nan')
        # case 2: There are valid values! Pick the first one in the array
        else:
            in_peaks[ column_index] = refl_col[ -1]


    # working on outbound data
    # since x_in and x_out are different values now :(
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( x_out), axis=0)

    out_refl = outbound_data.REFLECTIVITY.isel(time=0).isel(heading=0) # .transpose()

    # only for lat lon cases ( very annoying :/)
    if xaxis == 'lat' or xaxis == 'lon':
        # get rid of nan values from end of columns (like in plot_tdr() script in make_plots.py)
        out_refl = out_refl[ range( len( x_out) ), :]

    # turn all nans to zeros
    out_refl = np.where( np.isnan( out_refl.values), 0, out_refl)

    # return *indices* where our reflectivity value is greater than 20. other values are 0
    out_refl_ind = np.where( out_refl <= cutoff, 0, H_index_matrix )

    # make an empty array to hold eyewall location values
    out_peaks = np.empty( len( x_out), dtype=int)
    # cycle through each column until the first reflectivity value within cutoff is found: this is the start of the eyewall!
    for column_index in range( len( x_out)):
        refl_col = out_refl_ind[ column_index, :]

        # remove zeros: leave only indices with valid reflectivity returns
        refl_col = refl_col[ refl_col > 0]

        # case 1: no dBz values are high enough to meet threshold
        if np.size( refl_col) == 0:
            out_peaks[ column_index] = 0 # float('nan')
        # case 2: There are valid values! Pick the first one in the array
        else:
            out_peaks[ column_index] = refl_col[ -1]


    # renaming things
    H_in = H[in_peaks]
    H_out = H[out_peaks]

    # this code removes the zero height values found in the low dBz eye
    H_in = np.where( H_in > 0.0, H_in, float('nan') )
    H_out = np.where( H_out > 0.0, H_out, float('nan') )


    # optional code to cut off eyewall slope plotting at the first local maximum value!
    # the code defaults to the cutoff
    if eyewall_cutoff:
        in_local_max_vals = find_peaks( H_in, height = 8) # use the find peaks function to find local max values
        out_local_max_vals = find_peaks( H_out, height = 8) # use the find peaks function to find local max values

        if np.size( in_local_max_vals[0]) != 0:
            in_peak1_index = in_local_max_vals[0][0] # find the index of the first peak. the first [0] looks at all the indices, and the second [0] selects the first index
        else:
            in_peak1_index = 0
        if np.size( out_local_max_vals[0]) != 0:
            out_peak1_index = out_local_max_vals[0][0] # find the index of the first peak. the first [0] looks at all the indices, and the second [0] selects the first index
        else:
            out_peak1_index = 0

        x_in = x_in[ 0: in_peak1_index + 1] # the + 1 is to include the local peak
        H_in = H_in[ 0: in_peak1_index + 1]

        x_out = x_out[ 0: out_peak1_index + 1] # the + 1 is to include the local peak
        H_out = H_out[ 0: out_peak1_index + 1]

    return x_in, x_out, H_in, H_out




# another algorithm for finding the start of the eyewall, using in situ wind speed data! perhaps this is better?

# this needs a lot more work :/
def eyewall_start_in_situ(crl_path, crl_name, in_situ_path, in_situ_name, cutoff_indices, xaxis='time'):

    warnings.filterwarnings("ignore")

    # load data
    os.chdir( in_situ_path)
    xr_in_situ = plot_in_situ.load_in_situ( in_situ_path, in_situ_name, sample_step_size= 1)

    # rename variables from xarray for convenience
    float_time = xr_in_situ.float_time
    ws = [ float( line) for line in xr_in_situ["WS.d"].values] # convert strings to floats
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]

    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # deal with annoying indices :(
    if len( cutoff_indices) == 2:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
        time1 = crl_data.time[ i1]
        time2 = crl_data.time[ i2]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()
        # use the indices for nearest time values to trim down lat and lon values
        # this prevents data overlap and makes plotting faster!
        lon = lon[ idx1.values : idx2.values]
        lat = lat[ idx1.values : idx2.values]
        ws = ws[ idx1.values : idx2.values]

    elif len( cutoff_indices) == 4:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
        i3 = cutoff_indices[2]
        i4 = cutoff_indices[3]
        time1 = crl_data.time[ i1]
        time2 = crl_data.time[ i2]
        time3 = crl_data.time[ i3]
        time4 = crl_data.time[ i4]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()
        idx3 = (np.abs(float_time - time3)).argmin()
        idx4 = (np.abs(float_time - time4)).argmin()

        # use the indices for nearest time values to trim down lat and lon values
        # this prevents data overlap and makes plotting faster!
        lon = lon[ idx1.values : idx2.values] + lon[ idx3.values : idx4.values]
        lat = lat[ idx1.values : idx2.values] + lat[ idx3.values : idx4.values]
        ws = ws[ idx1.values : idx2.values] + ws[ idx3.values : idx4.values]
    else:
        print( 'Error in number of indices! Update them in tc_metadata.py')

    # pick the x axis to be plotted
    if xaxis == "time":
        xaxis_data = float_time
    elif xaxis == "lon":
        xaxis_data = lon
    elif xaxis == "lat":
        xaxis_data = lat
    else:
        print( "Please choose 'lon', 'lat', or 'time' as a valid xaxis!")


    # make sure the x axis data is a numpy array for proper plotting
    xaxis_data = np.array(xaxis_data)

    # find wind speed max points
    peaks = find_peaks( ws, prominence= 3, width=5, height= 15) # wlen=100

    # find minima to figure out where the center of the tc eye is
    # flip the wind speed dataset to use find peaks to get minima
    ws_flipped = [ -1 * val for val in ws]
    troughs = find_peaks( ws_flipped, prominence= 3, width=5, height= - 25)

    # find th, the height of the lowest trough (center of the eye)
    # th = trough height
    # - sign is to undo - sign above
    th = np.min( - troughs[1]['peak_heights'])
    # find the x position of the trough height
    # use xaxis_data[ ...] to go from an index to an x position
    # ti = trough index, tx = trough position
    ti = np.where( - troughs[1]['peak_heights'] == th)
    tx = xaxis_data[ troughs[0][ti]]
    # turn tx from a one element list into just that element
    tx = tx[0]

    # find the nearest left and right peaks to the trough tx!
    # left case:
    # find the peaks to the left of the tc eye center
    l_peaks = xaxis_data[ peaks [0]] [np.where( xaxis_data[peaks[0]] < tx)]
    # full list case (there are values to the left of the eye)
    if np.any( l_peaks):
        # lp(i, x, h) = left peak (index, position, height)
        lpi = (np.abs( l_peaks - tx )).argmin()
        lpx = l_peaks[ lpi]
        lph = ws[ np.where( xaxis_data == lpx) [0][0] ]
    # empty list case
    else:
        # just guess where the eyewall is lol... this needs to be fixed
        lpx = xaxis_data[ np.where( xaxis_data == tx)[0][0] - 50]
        lph = ws[ np.where( xaxis_data == lpx) [0][0] ]

    # right case
    r_peaks = xaxis_data[ peaks [0]] [np.where( xaxis_data[peaks[0]] >= tx)]
    # full list case
    if np.any( r_peaks):
        rpi = (np.abs( r_peaks - tx )).argmin()
        rpx = r_peaks[ rpi]
        rph = ws[ np.where( xaxis_data == rpx) [0][0] ]
    # empty list case
    else:
        rpx = xaxis_data[ np.where( xaxis_data == tx)[0][0] - 50]
        rph = ws[ np.where( xaxis_data == rpx) [0][0] ]

    in_x_start = lpx
    out_x_start = rpx

    warnings.filterwarnings("default")

    return in_x_start, out_x_start






# algorithm for finding the start of the eyewall using tdr data!
def eyewall_start( tdr_path, inbound_name, outbound_name, xaxis='dist'):

    warnings.filterwarnings("ignore")

    # using a cutoff of 20 dBz because it seems like other researchers have used this value for eyewall slope measurements!
    cutoff = 28.0 # i bumped this value up to 28 to get rid of effects from scud clouds in the eye
    # minimum height of a dBz point required to be considered the start of the eyewall (km)
    height_cutoff = 2.5 # height of P3 flight is around 3.5 km... would higher or lower work better?

    # load data
    os.chdir( tdr_path)
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # set up useful variables
    # choose a valid x axis
    if xaxis == 'dist':
        x_in = inbound_data.radius
        x_out = - outbound_data.radius
    elif xaxis == 'lon':
        # code taken from make_plots.py file, plot_tdr() function
        x_out = outbound_data.longitude [ ~np.isnan( outbound_data.longitude)]
        x_in =  inbound_data.longitude [ ~np.isnan( inbound_data.longitude)]
    elif xaxis == 'lat':
        # code taken from make_plots.py file, plot_tdr() function
        x_out = outbound_data.latitude  [ ~np.isnan( outbound_data.latitude)]
        x_in =  inbound_data.latitude  [ ~np.isnan( inbound_data.latitude)]
    else:
        print( 'choose \'lat\', \'lon\', or \'dist\' for an x axis! ')

    H = inbound_data.height
    H_index = range( len( H))
    # a matrix of height indices, to be used later. number of height values times number of x values
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( x_in), axis=0)


    # work on inbound data first
    in_refl = inbound_data.REFLECTIVITY.isel(time=0).isel(heading=0)
    # only for lat lon cases ( very annoying :/)
    if xaxis == 'lat' or xaxis == 'lon':
        # get rid of nan values from end of columns (like in plot_tdr() script in make_plots.py)
        in_refl = in_refl[ range( len( x_in) ), :]
    # turn all nans to zeros
    in_refl = np.where( np.isnan( in_refl.values), 0, in_refl)

    # return *indices* where our reflectivity value is greater than 20. other values are 0
    in_refl_ind = np.where( in_refl <= cutoff, 0, H_index_matrix )

    # cycle through each column until the first reflectivity value within cutoff is found: this is the start of the eyewall!
    for column_index in range( len( x_in)):
        refl_col = in_refl_ind[ column_index, :]

        # remove zeros: leave only indices with valid reflectivity returns
        refl_col = refl_col[ refl_col > 0]

        # make sure the current data meets the dBz and height requirements!
        if np.size( refl_col) != 0 and H[ refl_col[-1] ] > height_cutoff:
            # return the first valid element from the list
            in_x_start = x_in[ column_index]
            in_H_start = H[ refl_col[ -1]]
            break
        else:
            in_x_start = float('nan')
            in_H_start = float('nan')

    # working on outbound data
    H_index_matrix = np.repeat(np.array( H_index)[None, :], len( x_out), axis=0)
    out_refl = outbound_data.REFLECTIVITY.isel(time=0).isel(heading=0) # .transpose()
    # only for lat lon cases ( very annoying :/)
    if xaxis == 'lat' or xaxis == 'lon':
        # get rid of nan values from end of columns (like in plot_tdr() script in make_plots.py)
        out_refl = out_refl[ range( len( x_out) ), :]
    # turn all nans to zeros
    out_refl = np.where( np.isnan( out_refl.values), 0, out_refl)

    # return *indices* where our reflectivity value is greater than 20. other values are 0
    out_refl_ind = np.where( out_refl <= cutoff, 0, H_index_matrix )

    # cycle through each column until the first reflectivity value within cutoff is found: this is the start of the eyewall!
    for column_index in range( len( x_out)):
        refl_col = out_refl_ind[ column_index, :]
        # remove zeros: leave only indices with valid reflectivity returns
        refl_col = refl_col[ refl_col > 0]

        # make sure the current data meets the dBz and height requirements!
        if np.size( refl_col) != 0 and H[ refl_col[-1] ] > height_cutoff:
            # return the first valid element from the list
            out_x_start = x_out[ column_index]
            out_H_start = H[ refl_col[ -1]]
            break
        else:
            out_x_start = float('nan')
            out_H_start = float('nan')

    warnings.filterwarnings("default")

    return in_x_start, out_x_start, in_H_start, out_H_start
