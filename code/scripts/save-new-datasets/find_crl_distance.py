# import functions
import numpy as np
import os
import xarray as xr
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import helper_fns
import make_plots
os.chdir(  "/Users/etmu9498/research/code/scripts/save-new-datasets")
import save_crl_data
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts")
import load_in_situ_data
import clip_old_data




# this function finds a distance array for crl data using the tdr data!
def find_dist_new_tdr( tcname, dataset):
    # load data
    metadata = tc_metadata.all_data( tc= tcname)
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    tdr_path = "/Users/etmu9498/research/data/tdr-new"
    crl_path = "/Users/etmu9498/research/data/crl-new"
    os.chdir( tdr_path)
    new_tdr = xr.open_dataset( tdr_name)
    os.chdir( crl_path)

    new_crl = xr.open_dataset( crl_name)

    # goal 1: find the lat / lon limits for the crl dataset
    # using lat or lon?
    i1 = 0
    i2 = len( new_crl.time) - 1
    xtype = metadata[ 'xtype'][dataset]
    if xtype == 'lon':
        lim1 = new_crl.Lon[ i1]
        lim2 = new_crl.Lon[ i2]
        tdrx = new_tdr.longitude
    elif xtype == 'lat':
        lim1 = new_crl.Lat[ i1]
        lim2 = new_crl.Lat[ i2]
        tdrx = new_tdr.latitude
    else:
        print( 'update the xtype list in tc_metadata.all_data()!')

    # testing
    # print( xtype)
    # print( tdrx.values)
    # print( 'crl limit 1: ' + str( lim1.values))
    # print( 'crl limit 2: ' + str( lim2.values))

    # make a set of monatomically increasing values to act as an x axis
    inds = range( len( new_tdr.distance))

    # trying to fit slightly varying lat and lon tdr data with a polynomial!

    # get rid of nans because of course this is giving me trouble yet again...
    placeholder = -100000
    inds_nonan = np.where( ~ np.isnan( tdrx), inds, placeholder)
    inds_nonan = inds_nonan[ inds_nonan != placeholder]
    tdrx_nonan = tdrx[ ~ np.isnan( tdrx) ]

    # check output
    # print( inds)
    # print( tdrx)
    # print( len( inds))
    # print( len( tdrx))
    # trying to get rid of annoying nan values... this isn't working :/
    # print( len( tdrx.to_dataframe().dropna( inplace=True)))

    # fit the data without nans to a line!
    line_coefs = np.polyfit( inds_nonan, tdrx_nonan, 1)
    # print('line: y = ' + str( line_coefs[ 0]) + ' * x + ' + str( line_coefs[ 1]))

    # make a new line covering more values (aka extending the lat / lon range)
    # how many additional values should be added to the array?
    padding = 300
    # make a new index array to accomadate the extra values
    more_inds = np.linspace( 0 - padding, len( new_tdr.distance) + padding, num = len( new_tdr.distance) + 2*padding )
    # find the corresponding tdr lat / lon values for these extra distances!
    more_tdrx = more_inds * line_coefs[ 0] + line_coefs[ 1]

    # goal 2: try to find the closest lat / lon values in the tdr dataset
    # this is sooo much easier when only searching through one dataset!

    i1, x1 = helper_fns.closest_val( more_tdrx, lim1.values )
    i2, x2 = helper_fns.closest_val( more_tdrx, lim2.values )

    # print( 'tdr limit 1: ' + str( x1))
    # print( 'tdr limit 2: ' + str( x2))



    # goal 3: find the distances that correspond to the closest lat / lons found above!

    # repeat the process of making a larger dataset, but for distances this time
    inds_nonan = np.where( ~ np.isnan( new_tdr.distance), inds, placeholder)
    inds_nonan = inds_nonan[ inds_nonan != placeholder]
    dist_nonan = new_tdr.distance[ ~ np.isnan( new_tdr.distance) ]
    dist_coefs = np.polyfit( inds_nonan, dist_nonan, 1)

    # the new, extended distance array!
    more_dist = more_inds * dist_coefs[ 0] + dist_coefs[ 1]

    # print( more_dist[ i1-5:i1+5])
    # print( more_dist[ i2-5:i2+5])

    # setting the distance for the first crl limit
    dist1 = more_dist[ i1] # new_tdr.distance[ i1]
    # second limit
    dist2 = more_dist[ i2] # new_tdr.distance[ i2]

    # print( dist1)
    # print( dist2)

    '''
    testing!!!

    print( new_tdr.distance.values)
    print( 'tdr distance limit 1: ' + str( dist1))
    print( 'tdr distance limit 2: ' + str( dist2))

    print( 'length of distance array: ' + str( len( new_tdr.distance)))
    print( 'i1: ' + str( i1))
    print( 'i2: ' + str( i2))

    plt.scatter( inds, new_tdr.distance, s=5, c='g')

    plt.plot( more_inds, more_dist, c='y')
    plt.plot( range( len( new_tdr.distance) + 2 * padding), more_dist, c='b')
    plt.scatter( i1, dist1, s=50, c='b')
    plt.scatter( i2, dist2, s=50, c='b')

    plt.grid( 'on')
    plt.xlabel( 'Index')
    # plt.ylabel( 'Lat or Lon (Degrees)')
    plt.ylabel( 'Distance ( Km)')
    plt.show()
    plt.clf()
    '''

    # goal 4: make a distance array for the crl data given the distance tdr limits
    i1 = 0
    i2 = len( new_crl.time) - 1
    crl_distance = np.linspace( dist1, dist2, num= i2-i1)
    return crl_distance


# this function finds a distance array for crl data using in situ flight speed,
# time, and lowest p-3 flight height information
def find_crl_dist_in_situ( tcname, dataset, returnval='none'):

    warnings.filterwarnings("ignore")

    # load data
    metadata = tc_metadata.all_data( tc= tcname)

    in_situ_path = metadata['in_situ_path']
    in_situ_files = make_plots.load_flight_level( in_situ_path, print_files=False)
    in_situ_name = tc_metadata.choose_in_situ_date( metadata['dates'][dataset], in_situ_files)
    in_situ_data = load_in_situ_data.load_in_situ( in_situ_path, in_situ_name)


    if len( metadata['crl_range'][dataset]) == 4:
        crl_path = "/Users/etmu9498/research/data/crl-spiral-cases"
        crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])
    else:
        crl_path = "/Users/etmu9498/research/data/crl-new"
        tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)

    print( in_situ_name)
    print( crl_name)

    os.chdir( crl_path)
    new_crl = xr.open_dataset( crl_name)

    # time in hours
    float_time = in_situ_data.float_time
    # time in seconds (UTC)
    time = in_situ_data.float_time.values * 60 * 60
    # define the height of the P-3 above ocean surface
    height = [ float( line) for line in in_situ_data["HT.d"].values ]
    # true air speed of the P-3 (m/s)
    speed = [float( line) for line in in_situ_data["TAS.d"].values ]
    speed = np.array( speed)

    # interpolate between nans in speed array; they were causing issues!
    # old method
    # mask = np.isnan( speed)
    # speed[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), speed[~mask])
    nans, x = nan_helper( speed)
    speed[ nans] = np.interp( x( nans), x( ~nans), speed[~nans])

    # Part 1: Make a distance array for just the in situ dataset

    # This list will hold distance estimations for every in situ value, and it
    # will eventually be turned into a numpy array.
    # The first value of 0 establishes the baseline distance
    distance_array = [ 0.0]
    # do this for every data point in the numpy array, except the last value
    for index in range( len( height) - 1):

        speed_i = speed[ index]
        speed_i1 = speed[ index + 1]
        # most recent distance calculated
        dist_ind = distance_array[ -1]
        # average speed between indices
        avg_v = ( speed_i + speed_i1) / 2
        # distance gained during this time step
        deltax = avg_v * ( time[ index + 1] - time[ index])
        # make the new distance! and convert to km
        new_dist = dist_ind + deltax / 1000
        distance_array = distance_array + [ new_dist]

    # testing
    '''
    print( len( distance_array))
    print( np.nanmin( speed))
    print( np.nanmax( speed))
    print( speed[-20:-1])
    print( np.nanmin( time))
    print( np.nanmax( time))
    print( time[-20:-1])
    print( distance_array[-20:-1])
    '''

    # part two: making new distance axes with 0 centered on the lowest p-3 height!
    # this value corresponds to the lowest central pressure of the TC
    # total wind speed would also work well for this step

    # trim the data down to search the right section for the minimum value
    # maybe a little inefficient, but I know it'll work this way


    # trim height data down to the eye pass of interest
    trimheight = height_helper( new_crl, metadata['crl_range'][dataset], in_situ_data[ "HT.d"].values, float_time)

    '''
    print( distance_array[0:50])
    print( trimheight[0:50])
    '''
    
    # find the min height in this trimmed down dataset
    heightmin = np.nanmin( trimheight)
    heightmin_ind = np.nanargmin( trimheight) # this index is for the trimmed dataset

    # find the min height index for the full dataset!
    # the old way to find this index was flawed: it either just picked the first or
    # last occurence of the index and called that the min. but, this led to overlaps!
    # to fix this, we need to find the in situ index closest to the crl's starting index
    # and then add the trimheight index to it...
    idx1 = (np.abs(float_time.values - new_crl.time[ 0].values)).argmin()
    heightmin_fullind = heightmin_ind + idx1
    # old, wrong way
    # heightmin_fullind = np.where( height == heightmin)[0][-1]

    '''
    print( idx1)
    print( heightmin_ind)
    print('new index: ' + str( heightmin_fullind))
    print('old index ' + str( np.where( height == heightmin)[0][-1]))
    '''

    '''
    print( 'height:')
    print( heightmin)
    print( heightmin_ind)
    print( height[ heightmin_fullind])
    print( heightmin_fullind)
    '''

    # part two: find center distance value and subtract it to make it new 0 km
    center_dist = distance_array[ heightmin_fullind]
    # print( 'center distance: ' + str( center_dist))

    center_dist_array = np.subtract( np.array( distance_array), center_dist)

    '''
    print( distance_array[ heightmin_fullind - 10 : heightmin_fullind + 10])
    print( center_dist_array[ heightmin_fullind - 10 : heightmin_fullind + 10])
    '''

    # part 3: add new in situ distances to crl data
    time1 = new_crl.time[ 0]
    time2 = new_crl.time[ -1]

    # find the in situ times nearest the crl times
    idx1 = (np.abs(float_time - time1)).argmin()
    idx2 = (np.abs(float_time - time2)).argmin()


    # empty array that will hold new crl distances!
    crl_dist_array = []
    # new change 9/6/22
    # original case
    if len( metadata['crl_range'][dataset]) == 2:
        for crl_time in new_crl.time.values:
            # find the index of the in situ time nearest this crl time
            idx = (np.abs(float_time - crl_time)).argmin()
            # get the in situ distance value at this index
            crl_dist = center_dist_array[ idx]
            crl_dist_array = crl_dist_array + [ crl_dist]
    # new spiral case: do the same thing, but with the non chopped time array!
    elif len( metadata['crl_range'][dataset]) == 4:
        for crl_time in new_crl.time_distance_axis.values:
            idx = (np.abs(float_time - crl_time)).argmin()
            crl_dist = center_dist_array[ idx]
            crl_dist_array = crl_dist_array + [ crl_dist]

    '''
    # print distances and the times of the 0 km values for checking!
    print( len( crl_dist_array))
    print( np.min( crl_dist_array))
    print( np.max( crl_dist_array))
    print( 'In Situ Center')
    print( float_time[ heightmin_fullind].values)
    print( center_dist_array[ heightmin_fullind])
    print( float_time[ heightmin_fullind- 5: heightmin_fullind + 5].values)
    print( 'CRL center')
    idx = (np.abs( crl_dist_array)).argmin() # index of closest crl data point to 0 km
    print( new_crl.time[ idx].values)
    print( crl_dist_array[ idx])
    print( new_crl.time[ idx- 5: idx+5].values)
    # old code, this isn't necessary:
    # print( np.where( center_dist_array == 0.0))
    # Snp.all( [center_dist_array >= -0.0003, center_dist_array <= .0003 ] )))
    '''

    # center axis already accounts for 4 ind case! no need to change anything here
    if returnval == 'crl':
        return crl_dist_array

    # in situ data still doesn't account for the 4 ind case: let's fix that
    elif returnval == 'in-situ':
        # normal case: don't do anything
        if len( metadata['crl_range'][dataset]) == 2:
            return center_dist_array
        # 4 ind case: need to cut out distances during spiral
        if len( metadata['crl_range'][dataset]) == 4:
            # load trimmed crl data just for this bit to find spiral section
            temp_crl_path = "/Users/etmu9498/research/data/crl-new"
            os.chdir( temp_crl_path)
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
            temp_crl = xr.open_dataset( crl_name)
            # get time values
            newtime = temp_crl.time.values

            '''
            np.set_printoptions(threshold=np.inf)
            # difference between crl time axis with spiral removed and spiral kept!
            # for first 1/2 of array, time_diff = 0. When the spiral is removed, the
            # difference is around -.04 hours
            time_diff = new_crl.time_distance_axis.values - new_crl.time.values

            # first crl index with a time difference
            i1 = np.where( time_diff < - 0.001)[0][0]
            # last crl index with NO time difference
            i0 = i1 - 1
            '''
            # find inds where spiral starts / ends: crl data
            i0 = metadata['crl_range'][dataset][ 1] - metadata['crl_range'][dataset][ 0]
            i1 = metadata['crl_range'][dataset][ 2] - metadata['crl_range'][dataset][ 0]
            # find the indices of in situ times nearest the crl times
            # start and end in situ indices of spiral
            idx1 = (np.abs(float_time - newtime[i0])).argmin().values
            idx2 = (np.abs(float_time - newtime[i1])).argmin().values

            '''
            print(metadata['crl_range'][dataset])
            print( i0)
            print( i1)
            print( float_time[ idx1])
            print( newtime[i0])
            print( float_time[ idx2])
            print( newtime[i1])
            '''

            # shift end of array over through subtraction
            # figure out the distance difference between the two indices
            dist_diff = center_dist_array[idx1] - center_dist_array[idx2]

            # get the indices for the full distance array
            indrange = np.arange( len( center_dist_array))

            # for any index past the spiral, remove the distance difference to shift the axis over
            new_center_dist = np.where( indrange >= idx1, center_dist_array + dist_diff, center_dist_array)
            # remove in situ distances from the spiral region
            # new_center_dist = np.delete( new_center_dist, np.arange( idx1, idx2 )) # + 1))

            # add -999 values as fillers at the end of the dataset
            # new_center_dist = np.concatenate( ( new_center_dist, -999 * np.ones( idx2 - idx1)))

            warnings.filterwarnings("default")
            return new_center_dist

    else:
        print( "Please select either 'crl' or 'in-situ' as a valid input for" +
        " the returnval parameter in the find_crl_dist_in_situ() function.")
        return None



# this function finds a distance array for crl data using in situ flight speed,
# time, and lowest p-3 flight height information
# the center dist min value prevents a rmw from being selected right at 0 km! this creates
# a whole bunch of annoying problems
def find_crl_rmw( tcname, dataset, dist_lim=False, center_dist_min =10.0):

    warnings.filterwarnings("ignore")

    # load data
    metadata = tc_metadata.all_data( tc= tcname)

    # this in situ data should have spiral effects removed from it now!
    in_situ_path = metadata['new_flight_data_path']
    in_situ_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)
    os.chdir( in_situ_path)
    in_situ_data = xr.open_dataset( in_situ_name)

    # crl data with spirals removed- not the new height data, though!
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    os.chdir( metadata[ 'new_crl_path'])
    crl_data = xr.open_dataset( crl_name)


    # step 1: find in situ max winds for this tc eye

    # time in hours
    float_time = in_situ_data.float_time
    crltime = crl_data.time.values

    # trim in situ data down to crl range
    # in situ indices closest to crl time limits


    # dist_lim lets the user pick a distance range to look at the tc with!
    if dist_lim:
        # find in situ data points closest to dist lims: trim data to that size
        idx1 = (np.abs(in_situ_data.distance + dist_lim )).argmin().values
        idx2 = (np.abs(in_situ_data.distance - dist_lim )).argmin().values


        # make sure dist_lim fits within crl limits!
        # if not, just use the crl limits and print something out so the user knows this
        if in_situ_data.distance[ idx1] < crl_data.in_situ_distance[ 0] or in_situ_data.distance[ idx2] > crl_data.in_situ_distance[ -1]:
            idx1 = (np.abs(float_time - crltime[ 0])).argmin().values
            idx2 = (np.abs(float_time - crltime[ -1])).argmin().values
            print( "dist_lim value is too large: crl data limits will be used instead")

    # normal case: just use crl dataset limits
    else:
        idx1 = (np.abs(float_time - crltime[ 0])).argmin().values
        idx2 = (np.abs(float_time - crltime[ -1])).argmin().values

    # index where distance from center = 0 km in in situ data!
    # this used to be the old way of finding the center, but now a min distance is used
    idx3 = np.where( in_situ_data.distance == 0.0)[0][0]

    # find idx4 and idx5, the indices corresponding to -10 and 10 km in in situ data
    # indices for minimum distance to tc center (defaults to -10 and 10 km)
    idx4 = (np.abs(in_situ_data.distance + center_dist_min )).argmin().values # -10
    idx5 = (np.abs(in_situ_data.distance - center_dist_min )).argmin().values # 10 km

    # then, edit the code below to split up wind speeds into positive and negative
    # distances from -10 and 10 km, not 0 km!
    # or, use find_peaks() to do this? more control on which peaks are legit, etc?
    if idx1 < idx2:
        ws_negative = in_situ_data['WS.d'][idx1:idx4]
        ws_positive = in_situ_data['WS.d'][idx5+1:idx2]
    else:
        ws_negative = in_situ_data['WS.d'][idx2:idx4]
        ws_positive = in_situ_data['WS.d'][idx5+1:idx1]


    ''' original code: center = 0 km
    # make sure all wind speed values are floats
    # trimming down the wind speed dataset first saves a bunch of time!
    # break wind speed data down into before and after 0 km mark
    if idx1 < idx2:
        ws_negative = in_situ_data['WS.d'][idx1:idx3]
        ws_positive = in_situ_data['WS.d'][idx3+1:idx2]
    else:
        ws_negative = in_situ_data['WS.d'][idx2:idx3]
        ws_positive = in_situ_data['WS.d'][idx3+1:idx1]

    this code came after float conversion line... inds were changed below

    # find inds of max wind values- one for each side of 0 km
    # add idx1 to get back to original indices
    if idx1 < idx2:
        i_max1 = np.nanargmax( ws_negative) + idx1
        i_max2 = np.nanargmax( ws_positive) + idx3 + 1
    else:
        i_max1 = np.nanargmax( ws_negative) + idx2
        i_max2 = np.nanargmax( ws_positive) + idx3 + 1


    '''

    ws_negative = [ float( value) for value in ws_negative]
    ws_positive = [ float( value) for value in ws_positive]

    # find inds of max wind values- one for each side of 0 km
    # add idx1 to get back to original indices
    if idx1 < idx2:
        i_max1 = np.nanargmax( ws_negative) + idx1
        i_max2 = np.nanargmax( ws_positive) + idx5 + 1
    else:
        i_max1 = np.nanargmax( ws_negative) + idx2
        i_max2 = np.nanargmax( ws_positive) + idx5 + 1


    ''' testing
    this code highlights the issues in picking 0 km as the center and finding maxes from there!

    print( in_situ_data.distance[ i_mindist1].values)
    print( in_situ_data.distance[ i_mindist2].values)
    print( in_situ_data.distance[ i_max1].values)
    print( in_situ_data.distance[ i_max2].values)

    # make sure both dist_lim values are outside of 10 km!
    # if they aren't, repeat the process above but with max values outside of 10 km
    if in_situ_data.distance[ i_mindist1] < in_situ_data.distance[ i_max1] or in_situ_data.distance[ i_mindist2] > in_situ_data.distance[ i_max2]:
        idx1 = i_mindist1
        idx2 = i_mindist2
        print( "Max wind speed value is too close to the tc center: a limit of " + str(center_dist_min) + ' will be used instead')
    '''

    # convert in situ indices to crl indices
    # these indices hold rmw = 1 values!
    crl_idx1 = (np.abs(crltime - float_time[ i_max1].values)).argmin()
    crl_idx2 = (np.abs(crltime - float_time[ i_max2].values)).argmin()

    # speeds and distances look good! axes are still misaligned, maybe it's a
    # problem with the slope determinations?
    print( 'speed 1: ' + str( in_situ_data['WS.d'][ i_max1].values))
    print( 'speed 2: ' + str( in_situ_data['WS.d'][ i_max2].values))

    print( 'distance 1: ' + str( crl_data.in_situ_distance[ crl_idx1].values))
    print( 'distance 2: ' + str( crl_data.in_situ_distance[ crl_idx2].values))


    ''' testing
    print( in_situ_data.distance[ idx1].values)
    print( in_situ_data.distance[ idx2].values)

    print( idx1)
    print( idx3)
    print( idx2)

    print( in_situ_data['WS.d'][i_max1].values)
    print( in_situ_data['WS.d'][i_max2].values)
    print( in_situ_data['distance'][i_max1].values)
    print( in_situ_data['distance'][i_max2].values)

    print( crl_data.in_situ_distance[ crl_idx1].values)
    print( crl_data.in_situ_distance[ crl_idx2].values)
    print( in_situ_data['distance'][idx3].values)
    '''

    # define location of max in situ winds as rmw = 1

    # make a linear function to convert distance to rmw
    # y intercept = 0 bc rmw = 0 when distance = 0km
    # so, just find the slope of the line
    # account for negative distances; rmw is always positive
    if crl_data.in_situ_distance[ crl_idx1] < 0:
        slope1 = 1 / crl_data.in_situ_distance[ crl_idx1].values
    else:
        slope1 = - 1 / crl_data.in_situ_distance[ crl_idx1].values
    # repeat for ind2
    if crl_data.in_situ_distance[ crl_idx1] < 0:
        slope2 = 1 / crl_data.in_situ_distance[ crl_idx2].values
    else:
        slope2 = - 1 / crl_data.in_situ_distance[ crl_idx2].values

    # create and return an rmw axis!
    # split crl distances into positive and negative segments
    dist_pos = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance > 0.0)].values
    dist_neg = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance <= 0.0)].values

    # apply the slopes to find rmw values!
    rmw_pos = slope2 * dist_pos
    rmw_neg = slope1 * dist_neg

    rmw = np.concatenate( (rmw_neg, rmw_pos))

    # in addition to the totally positive rmw variable, return an rmw axis with negative values!
    # this will help with plotting crl data later
    slope1_neg = - 1 / crl_data.in_situ_distance[ crl_idx1].values
    slope2_neg =  1 / crl_data.in_situ_distance[ crl_idx2].values
    rmw_2 = slope2_neg * dist_pos
    rmw_1 = slope1_neg * dist_neg
    rmw_negative = np.concatenate( (rmw_1, rmw_2))


    print( rmw_negative[ crl_idx1])
    print( rmw_negative[ crl_idx2])

    '''
    # check for infs or nans
    count = 0
    count2 = 0
    for item in rmw:
        if np.isnan( item):
            count += 1
        if np.isinf( item):
            count2 += 1

    if count > 0 or count2 > 0:
        print( 'number of nans: ' + str( count))
        print( 'number of infs: ' + str( count2))
        print( 'slope 1: ' + str( slope1))
        print( 'slope 2: ' + str( slope2))
    '''

    '''
    # some nans seem to be sneaking through; interpolate to get rid of them
    nans, x = nan_helper( rmw)
    rmw[ nans] = np.interp( x( nans), x( ~nans), rmw[~nans])

    nans, x = nan_helper( rmw_negative)
    rmw_negative[ nans] = np.interp( x( nans), x( ~nans), rmw_negative[~nans])
    '''

    return rmw, rmw_negative




# helper function to replace nans in distance array! taken from stack overflow lol
def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]



# helper function to convert an in situ return variable into a float and
# find its corresponding crl data point
def height_helper(crl_data, cutoff_indices, return_var, float_time):
    return_var_temp = np.zeros( len( return_var))
    for line_ind in range( len( return_var)):
        if return_var[ line_ind] == '':
            return_var_temp[line_ind] = np.nan
        else:
            return_var_temp[ line_ind] = float( return_var[ line_ind])
    return_var = return_var_temp.tolist()

    # deal with annoying indices :(
    # normal case
    if len( cutoff_indices) == 2:
        time1 = crl_data.time[0]
        time2 = crl_data.time[-1]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()

    elif len( cutoff_indices) == 4:
        # use special time axis for these cases to avoid time jump!
        time1 = crl_data.time_distance_axis[0]
        time2 = crl_data.time_distance_axis[-1]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()

    else:
        print( 'error error error!')

    return return_var[ idx1.values : idx2.values]
