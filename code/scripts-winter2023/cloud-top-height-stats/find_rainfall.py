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


# this function uses find_peaks() to find rainfall peak values!
# This is used as a helper function in find_rain(), which finds and returns the
# actual full regions of rain
def find_rain_peaks( power, cutoff_power, xaxis, H_index_matrix, p3_heights, H, grace_case, cloud_heights=False):

    # filter the power and power indices from H_index_matrix
    power = power.where( power.values > cutoff_power)
    # this line changes all NaNs to 0s, and all non Nans to indices
    power_ind_matrix = np.where( np.isnan( power) , 0, H_index_matrix)

    # get rid of nans in the power matrix, they were causing errors for some reason.
    # -200 is far below every max peak.
    # not sure why I need to do this really
    power = np.where( np.isnan( power) , -200, power)

    # this empty list will be returned at the end of the function call
    # the list allows us to return multiple height values
    # using lists instead of numpy arrays here for efficiency
    power_index = []

    # cycle through each column and find cloud heights
    for column_index in range( len( xaxis)):
        # select the current column from power and power_ind_matrix. Then, get rid of all 0s
        power_column = power[ column_index, :]
        power_ind_column = power_ind_matrix[ column_index, :]
        power_ind_column = power_ind_column[ power_ind_column > 0]


        # check for empty lists for both power_inds and powers: no backscattered values at all
        # in this case, there should be no return height
        if np.size( power_ind_column) == 0:

            # new code: return an empty list
            power_index += [ []]
            '''
            # old way of doing things...
            # wrap all index values with parentheses to differentiate between the different runs!
            power_index += [ [0]]
            '''
            continue

        if np.size( power_column) == 0:
            power_index += [ []]
            continue


        # otherwise, only look at the first clear air patch
        # split current list of indices into arrays of consecutive values! Helps determine the top cloud chunk
        splits = np.split( power_ind_column, np.where(np.diff( power_ind_column) != 1)[0] + 1)
        # only look at the first clear air patch
        first_consec_inds = splits[0]
        # trim the power column to the size of the first clear air patch
        power_column = power_column[ first_consec_inds[0] : first_consec_inds[-1] ]

        # new code: remove values at and below cloud heights
        if cloud_heights:
            # the [0] is to look at the top cloud, if there are multiple clouds
            power_column = power_column[ 0 : cloud_heights[ column_index][0] ]


        # width makes sure that small, sharp peaks are not included
        #     -> this value is larger than cloud top find_peaks()
        # prominence makes sure that smooth peaks are significant relative to one another
        #     -> this value is smaller than cloud top find_peaks()
        # height makes sure that a minimum backscatter threshold is met
        #     -> without this, clear sky regions might be incorrectly classified as rain
        # I decided to get rid of the wlen parameter now, because it doesn't matter how
        # wide the upper limit of the rainfall signal is!

        '''
        print( 'len of power for index ' + str( column_index) + ": " + str( len( power_column)))
        print( power_column)
        '''

        peaks = find_peaks( power_column, prominence= 1, height=-17.5, width = 25 )

        # new code: if no peaks exist, just append an empty list to the return values!
        # no valid rain layers
        if np.size( peaks[0]) == 0:
            power_index += [ []]

            '''
            # old way of doing things...
            # if no peaks exist, pick the end of the first clear cloud layer
            if np.size( peaks[0]) == 0:

                # annoying grace cases
                # these lines might actually be wrong! maybe do power_index += instead
                if grace_case == 1:
                    ind_first_power_val = ( np.abs( np.subtract( H, p3_heights[ column_index] + .4))).argmin()
                    power_index[ column_index] = first_consec_inds[ -1] - ind_first_power_val
                elif grace_case == 2:
                    ind_first_power_val = ( np.abs( np.subtract( H, p3_heights[ column_index] - .65))).argmin()
                    power_index[ column_index] = first_consec_inds[ -1] - ind_first_power_val
                # normal case
                else:
                    # look at the last value in the first array of non NaN values: this represents
                    # the end of the first clear air chunk beneath the plane! Aka the cloud top
                    # this value actually needs to get shifted up by the distance of the blank
                    # space! This line of code was causing some problems with the new in situ data
                    ind_first_power_val = ( np.abs( np.subtract( H, p3_heights[ column_index]))).argmin()

                    # shift values up by the distance contained in the blank space!
                    # power_index[ column_index] = first_consec_inds[ -1] - ind_first_power_val
                    power_index += [[ first_consec_inds[ -1] - ind_first_power_val.values]]
            '''

        # case where local peaks exist! Use the prominence output from find_peaks
        # one or more peaks exist: return them all!
        else:
            # peaks[0] returns all peaks, peaks[0][0] would return only the first peak
            power_index += [ peaks[0] ]

    # print( power_index)
    return power_index


# This wrapper function uses find_rain_peaks() and some more data manipulation to
# identify regions of rainfall
def find_rain( power, cutoff_power, xaxis, H_index_matrix, p3_heights, H, grace_case, cloud_heights=False, rain_cutoff=-20):

    power_index = find_rain_peaks( power, cutoff_power, xaxis, H_index_matrix, p3_heights, H, grace_case, cloud_heights)


    # filter the power and power indices from H_index_matrix
    power = power.where( power.values > rain_cutoff)

    '''
    # this line changes all NaNs to 0s, and all non Nans to indices
    power_ind_matrix = np.where( np.isnan( power) , 0, H_index_matrix)

    # get rid of nans in the power matrix, they were causing errors for some reason.
    # -200 is far below every max peak.
    # not sure why I need to do this really
    power = np.where( np.isnan( power) , -200, power)
    '''

    # this empty list will be returned at the end of the function call
    # the list allows us to return multiple height values
    # using lists instead of numpy arrays here for efficiency
    rain_index = []

    # cycle through each column and find cloud heights
    for column_index in range( len( xaxis)):
        # select the current column from power and power_ind_matrix. Then, get rid of all 0s
        power_column = power[ column_index, :]

        '''
        power_ind_column = power_ind_matrix[ column_index, :]
        power_ind_column = power_ind_column[ power_ind_column > 0]
        '''

        # maybe add splits code here if things aren't working?

        # remove values at and below cloud heights
        if cloud_heights:
            # the [0] is to look at the top cloud, if there are multiple clouds
            power_column = power_column[ 0 : cloud_heights[ column_index][0] ]

        '''
        # no rainfall cases:empty power inds
        if np.size( power_ind_column) == 0:
            rain_index += [ []]
            continue
        '''

        # empty powers
        if np.size( power_column) == 0:
            rain_index += [ []]
            continue

        # no rainfall peak
        if np.size( power_index[ column_index]) == 0:
            rain_index += [ []]
            continue
        # add rainfall peak here!
        else:

            # turn power values into indices and remove nans
            power_column = np.where( np.logical_not( np.isnan( power_column) ))[0].tolist()
            rain_index += [ power_column ]

    # print( '\n\n\n')
    # print( rain_index)
    return rain_index
