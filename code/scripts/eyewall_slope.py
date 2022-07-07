import os
import xarray as xr
import numpy as np
from scipy.signal import find_peaks

os.chdir( "/Users/etmu9498/research/code/scripts/")
import make_plots



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
        '''
        print( 'inbound size with no nans: ' + str( np.size( inbound_data.longitude[ ~np.isnan( inbound_data.longitude)], 0)) )
                # + ' x ' + str( np.size( inbound_data.longitude[ ~np.isnan( inbound_data.longitude)], 1)) )
        print( 'inbound size with nans: ' + str( np.size( inbound_data.longitude, 0)) ) # + ' x ' + str( np.size( inbound_data.longitude, 1)) )

        print( 'outbound size with no nans: ' + str( np.size( outbound_data.longitude[ ~np.isnan( outbound_data.longitude)], 0)))
                # + ' x ' + str( np.size( outbound_data.longitude[ ~np.isnan( outbound_data.longitude)], 1)))
        print( 'outbound size with nans: ' + str( np.size( outbound_data.longitude, 0)) ) # + ' x ' + str(np.size( outbound_data.longitude, 1)))
        '''

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





# algorithm for finding points along the eyewall!
def eyewall_start( tdr_path, inbound_name, outbound_name, xaxis='dist'):

    # using a cutoff of 20 dBz because it seems like other researchers have used this value for eyewall slope measurements!
    cutoff = 20.0
    # minimum height of a dBz point required to be considered the start of the eyewall (km)
    height_cutoff = 3.5 # height of P3 flight... would higher or lower work better?

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


    return in_x_start, out_x_start, in_H_start, out_H_start
