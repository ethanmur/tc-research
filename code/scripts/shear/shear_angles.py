import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

os.chdir( "/Users/etmu9498/research/code/scripts")
import tc_metadata


# using the environmental shear data saved in the tc_metadata.py file and the determined
# p3 flight direction, figure out which shear quadrants are crossed. The first returned
# value is the quadrant entered, while the second is the exit quadrant!
def two_shear_quadrants( tcname, dataset, close_warn=False):
    # load metadata and dataset names for this case
    metadata = tc_metadata.all_data( tcname)
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    in_situ_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

    print( in_situ_name)

    # load data path info
    crl_path = metadata[ 'new_crl_path']
    in_situ_path = metadata[ 'new_flight_data_path']

    # load data
    os.chdir( in_situ_path)
    in_situ_data = xr.open_dataset( in_situ_name)

    # get the environmental shear for this example
    shear_dir = metadata['shear_dir'][ dataset]

    # find the p-3 heading at the TC center without the helper script
    # get the in situ distance from tc center axis
    dist = in_situ_data.distance
    # find the center (distance = 0 km)
    zero_ind = np.where( dist == 0.0)
    # find the flight header at the tc center
    p3_dir = float( in_situ_data.THDGref[ zero_ind].values[ 0])

    # figure out the quadrants! do this through angles and the p-3 position
    difference = p3_dir - shear_dir

    if ( difference >= 0 and difference < 90) or ( difference >= -360 and difference < -270):
        early_quad = 'UL'
        late_quad = 'DR'
    elif (difference >= 90 and difference < 180) or ( difference >= -270 and difference < -180):
        early_quad = 'DL'
        late_quad = 'UR'
    elif (difference >= 180 and difference < 270) or ( difference >= -180 and difference < -90):
        early_quad = 'DR'
        late_quad = 'UL'
    elif (difference >= 270 and difference < 360) or ( difference >= -90 and difference < 0):
        early_quad = 'UR'
        late_quad = 'DL'
    else:
        print( 'error! something is wrong with the if statement in two_shear_quadrants()')

    # Optional code to test to see if the flight path is within 10 degrees of a shear quadrant limit
    # If so, print a warning for the user
    if close_warn:

        abs_angle = np.abs( difference)
        print( abs_angle)

        quad_angles = [ 0, 90, 180, 270]
        for qa in quad_angles:
            if abs_angle - qa < 10 and abs_angle - qa > 0:
                print( "Warning: Flight path is within 10 degrees of shear quadrant boundary.")

    return early_quad, late_quad


def in_situ_p3_dir( tcname, dataset):
    # load metadata and dataset names for this case
    metadata = tc_metadata.all_data( tcname)
    in_situ_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

    # load data path info
    in_situ_path = metadata[ 'new_flight_data_path']
    # load data
    os.chdir( in_situ_path)
    in_situ_data = xr.open_dataset( in_situ_name)

    # get the in situ distance from tc center axis
    dist = in_situ_data.distance
    # find the center (distance = 0 km)
    zero_ind = np.where( dist == 0.0)
    # find the flight header at the tc center
    p3_angle = float( in_situ_data.THDGref[ zero_ind].values[ 0])
    return p3_angle

    '''
    print( 'tc ' + tcname + ' dataset ' + str( dataset))
    print( 'distance at tc center: ' + str( dist[ zero_ind].values))
    print( "P-3 heading at TC center: " + str( in_situ_data.THDGref[ zero_ind].values) + '\n')

    print( 'number of distance points: ' + str( len( dist)))
    print( 'number of distance = 0 km points: ' + str( len( zero_ind)))
    print( dist[ zero_ind])
    print( zero_ind)
    '''
def metadata_shear_dir( tcname, dataset):
    # load metadata and dataset names for this case
    metadata = tc_metadata.all_data( tcname)
    return metadata['shear_dir'][dataset]


# this function finds the angle of the P-3 flight path. It uses wind shear conventions
# when finding angles: 0 degrees points from the South, 45 degrees from the SE, etc...
def plot_flight_angles( tcname, dataset):
    # load metadata and dataset names for this case
    metadata = tc_metadata.all_data( tcname)
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    in_situ_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

    # load data path info
    crl_path = metadata[ 'new_crl_path']
    in_situ_path = metadata[ 'new_flight_data_path']

    # load data
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # clip lat and lon data to only include the straight line flight track segment
    shear_ind0 = metadata[ 'crl_shear_inds'][ dataset][ 0]
    shear_ind1 = metadata[ 'crl_shear_inds'][ dataset][ 1]
    data_len = len( crl_data.time)
    lat = crl_data.Lat[ 0 + shear_ind0: data_len - shear_ind1]
    lon = crl_data.Lon[ 0 + shear_ind0: data_len - shear_ind1]

    # fit the lat lon flight line to eveutally figure out
    line_coefs = np.polyfit( lon, lat, 1)

    '''
    print( 'time length (no clipping): ' + str( len( crl_data.time)))
    print( 'ind 1: ' + str( shear_ind0))
    print( 'ind 2: ' + str( shear_ind1))
    print( 'lat length: ' + str( len( lat)))
    print( 'lon length: ' + str( len( lon)))
    print('line: y = ' + str( line_coefs[ 0]) + ' * x + ' + str( line_coefs[ 1]))
    '''

    # make a new line covering more values (aka extending the lat / lon range)
    # how many additional values should be added to the array?
    low_lim = np.nanmin( lon) - 1 # -62
    high_lim = np.nanmax( lon) + 1 # -48
    x = np.linspace( low_lim, high_lim, num=100)
    y = line_coefs[ 0] * x + line_coefs[ 1]

    # plot results!
    '''
    plt.figure( figsize=(12, 6))
    plt.scatter( crl_data.Lon, crl_data.Lat, s=3, c='g')
    plt.plot( x, y, c='y')
    plt.scatter( lon[ 0], lat[0], s=15, c='b')
    plt.scatter( lon[ -1], lat[ -1], s=15, c='b')

    plt.grid( 'on')
    plt.xlabel( 'Lon')
    plt.ylabel( 'Lat')

    plt.xlim( [ low_lim, high_lim])
    plt.ylim( [ np.nanmin( lat) - 1, np.nanmax( lat) + 1])
    plt.show()
    '''
