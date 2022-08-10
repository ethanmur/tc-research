# import functions
import numpy as np
import os
import xarray as xr

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
    # load data
    metadata = tc_metadata.all_data( tc= tcname)

    in_situ_path = metadata['in_situ_path']
    in_situ_files = make_plots.load_flight_level( in_situ_path, print_files=False)
    in_situ_name = tc_metadata.choose_in_situ_date( metadata['dates'][dataset], in_situ_files)
    in_situ_data = load_in_situ_data.load_in_situ( in_situ_path, in_situ_name)

    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    crl_path = "/Users/etmu9498/research/data/crl-new"
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

    '''
    print( len( time))
    print( len( height))
    print( len( speed))
    '''

    # Part 1: Make a distance array for just the in situ dataset

    # This list will hold distance estimations for every in situ value, and it
    # will eventually be turned into a numpy array.
    # The first value of 0 establishes the baseline distance
    distance_array = [ 0.0]
    # do this for every data point in the numpy array, except the last value
    for index in range( len( height) - 1):
        # most recent distance calculated
        dist_ind = distance_array[ -1]
        # average speed between indices
        avg_v = ( speed[ index] + speed[ index + 1]) / 2
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

    # trim height data
    crl_old_names = make_plots.load_crl( metadata['crl_path'], print_files=False)
    crl_old_name = tc_metadata.choose_crl_date( metadata['dates'][dataset], crl_old_names)
    cutoff_indices = metadata['crl_range'][dataset]
    trimheight = clip_old_data.in_situ_helper( metadata['crl_path'], crl_old_name, cutoff_indices, in_situ_data[ "HT.d"].values, float_time)

    # find the min height in this dataset
    heightmin = np.nanmin( trimheight)
    heightmin_ind = np.nanargmin( trimheight) # this index is for the trimmed dataset
    heightmin_fullind = np.where( height == heightmin)[0][0] # this index is for the full dataset!

    '''
    print( 'height:')
    print( heightmin)
    print( heightmin_ind)
    print( height[ heightmin_fullind])
    print( heightmin_fullind)
    '''

    # part two: find center distance value and subtract it to make it new 0 km
    center_dist = distance_array[ heightmin_fullind]
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
    for crl_time in new_crl.time.values:
        # find the index of the in situ time nearest this crl time
        idx = (np.abs(float_time - crl_time)).argmin()

        # get the in situ distance value at this index
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

    if returnval == 'crl':
        return crl_dist_array
    elif returnval == 'in-situ':
        return center_dist_array
    else:
        print( "Please select either 'crl' or 'in-situ' as a valid input for" +
        " the returnval parameter in the find_crl_dist_in_situ() function.")
        return None
