# import functions
import numpy as np
import os
import xarray as xr
import warnings
import scipy

os.chdir(  "/Users/etmu9498/research/code/scripts")
import helper_fns
import make_plots




# this function finds a matching, processed flight level data file for transferring
# rmw, dists, and p-3 height calculations, among other useful fields!
def find_matching_fl( year, crlname):
    # first, go to the correct year folder for the crl / fl case
    flpath = "/Users/etmu9498/research/data/in-situ-noaa-processed/" + year

    # get a list of the files within this folder
    fllist = make_plots.load_flight_level( flpath, print_files=False)

    # go through all of the files within fllist. Save all the ones with matching dates
    # with the crl data in a new list!
    matchlist = []
    for fl_filei, fl_fileval in enumerate( fllist):
        if crlname[3:-18] == fl_fileval[ : 10]:
            print( fl_fileval[ : 10])
            matchlist.append( fl_fileval)

    # error cases: either too many files with the same date from the same plane, or no valid files!
    if len( matchlist) == 0:
        print(crlname[3:-18])
        print( "Error! No flight level dataset found. Investigate!")
        return False
    elif len( matchlist) > 1:
        print( "Error! Too many flight level datasets found for the same date and aircraft. Investigate!")
        print(crlname[3:-18])
        return False

    # get the fl dataset and return it!
    os.chdir( flpath)
    return xr.open_dataset( matchlist[ 0], decode_times=True)




# downscale the existing flight level distance and rmw axes and return them for saving
# in the new crl dataset!
def find_rmws( crl_data, fl_data):
    # get the flight level distance and rmw axes
    fldist, flrmw = fl_data['center_dist'], fl_data['rmw']
    #print("first crl time: " + str( crl_data.time[0].values))
    #print("last crl time: " + str( crl_data.time[-1].values))
    #print("first fl time: " + str( fl_data['float_time'][0].values))
    #print("last fl time: " + str( fl_data['float_time'][-1].values))

    # trim down the distances and rmws to the crl's time limits!
    # flight level indices of min + max crl times!
    i0 = np.argmin( np.abs( fl_data.time.values - crl_data.time[0].values ))
    i1 = np.argmin( np.abs( fl_data.time.values - crl_data.time[-1].values ))

    # get the flight level distance and rmw axes
    fldist, flrmw = fldist[ i0:i1], flrmw[i0:i1]

    # interpolate these values down to the crl's matrix size!
    # older numpy method:
    # this method results in nearly identical output as the scipy function below!
    # values are sligtly shifted here for some reason
    # feel free to switch to this method if scipy is misbehaving
    # crldist = np.interp( crl_data['time'], fl_data['float_time'][i0:i1], fldist)

    dist_interp = scipy.interpolate.interp1d( np.arange(fldist.size), fldist)
    crldist = dist_interp( np.linspace( 0, fldist.size - 1, crl_data.time.size))

    rmw_interp = scipy.interpolate.interp1d( np.arange( flrmw.size), flrmw)
    crlrmw = rmw_interp( np.linspace( 0, flrmw.size - 1, crl_data.time.size))
    return crldist, crlrmw


# downscale the existing chosen flight level value onto a crl scale! This works a lot
# like find_rmws() but obviously focusing on correcting heights. Really useful when
# correcting / interpolating crl matrices
def find_interp( crl_data, fl_data, fl_data_field):

    # trim down the distances and rmws to the crl's time limits!
    # flight level indices of min + max crl times!
    i0 = np.argmin( np.abs( fl_data.time.values - crl_data.time[0].values ))
    i1 = np.argmin( np.abs( fl_data.time.values - crl_data.time[-1].values ))
    fl_data_field = fl_data_field[ i0:i1]

    interp = scipy.interpolate.interp1d( np.arange(fl_data_field.size), fl_data_field)
    crl_field = interp( np.linspace( 0, fl_data_field.size - 1, crl_data.time.size))
    return crl_field


# interpolate crl matrices using the updated p3 height measurements!
# based on code from scripts/in-situ-scripts/get_p3_heights.py (interp_data2)
def interp_data( varval, crl_origh, p3_heights):

    resolution = 6 # change this?

    # defining a new base height array to standardize the new matrix
    base_h = np.arange( -4.0, 0.0, resolution / 1000)
    # make an empty array for initial storage
    new_matrix=np.empty([ np.size( varval, 0), len( base_h)])


    # redefine variable
    origh = crl_origh
    orig_maxh = np.nanmax( - origh)

    # do this for every x axis value
    for i in range( np.size( varval, 0) ):
        # get the current matrix column
        columni = varval[ i, :]

        # calculate a couple more useful heights: height from p3 and difference
        # from original estimate
        p3h = p3_heights[ i] / 1000
        hdiff = orig_maxh - p3h

        # make a new height profile for this column, from 0 to p-3 height
        # this is different from the step above!
        # ex: - p3h ~ -3.15, 0 - hdiff ~ + .4, should get rid of lower vals!
        trim_h = np.linspace( - p3h , 0.0 + hdiff, np.size( varval, 1))


        # interpolate over that row and save results
        new_columni = np.interp( base_h, trim_h, columni)
        new_matrix[ i, :] = new_columni

    return base_h, new_matrix


def find_rh_old( crl_data, fl_data):
    # relative humidity
    # defining constants
    epsilon = .622
    e_0 = 6.112 # hPa
    b = 17.67
    T_1 = 273.15 # constants (in K)
    T_2 = 29.65 # K
    temp = xr.DataArray( T_2d.transpose()) + 273 # temperature in kelvin
    wvmr = xr.DataArray( wv_2d.transpose())

    # define pressure field using scale height
    # no negative sign in exp() because heights are already negative
    scale_ht = 7.5 # km, just an estimate
    pressure = 1013.3 * np.exp( xr.DataArray( newh) / scale_ht) # hPa

    # find saturation vapor pressure
    e_s = e_0 * np.exp(  ( b * ( temp - T_1) ) / (temp - T_2) )

    # find relative humidity!
    saturation_wvmr = 1000 * (epsilon * e_s ) / ( pressure - e_s) # multiply by 1000 to get to g/kg like wvmr (above)
    crl_rh = 100 * wvmr / saturation_wvmr

    '''
    print( 'temp: ' + str( np.shape( temp)))
    print( 'wvmr: ' + str( np.shape( wvmr)))
    print( 'pressure: ' + str( np.shape( pressure)))
    print( 'es: ' + str( np.shape( e_s)))
    print( 'sat mr: ' + str( type( saturation_wvmr)))
    print( 'sat mr: ' + str( np.shape( saturation_wvmr)))
    print( 'rh: ' + str( type( crl_rh)))
    print( 'rh: ' + str( np.shape( crl_rh)))
    '''

    crl_rh = crl_rh.where( crl_rh.values <= 100)
    crl_rh = crl_rh.where( crl_rh.values >= 0)

    # try going from xr to list back to xr?
    # transpose flips things back to correct order
    # this is all to avoid an error while saving the matrix later... and it
    # seems to work!
    crl_rh = xr.DataArray( crl_rh).transpose()
    crl_rh = crl_rh.to_numpy()
    crl_rh = crl_rh.tolist()
    crl_rh = np.array( crl_rh)
    crl_rh = xr.DataArray( crl_rh)


def find_temp_anom_old( crl_data, fl_data):
    # temperature anomaly
    # temperature anomaly is found using a local average from our data! An array
    layer_avg_temp = np.nanmean( T_2d, axis = 0)
    # make an array of average temperatures within the correct bounds. A matrix
    avg_temp_obs = np.ones_like( T_2d)


    # for every height value:
    for i in range( np.size( T_2d, 0)):
        # since avg_temp_obs is initially all 1s, multiplying by layer_avg_temp
        # just gives it that value quickly!
        avg_temp_obs[i, :] = avg_temp_obs[ i, :] * layer_avg_temp
    # calculate the anomaly: current value - layer average for every layer
    temp_anomaly_obs = T_2d - avg_temp_obs
