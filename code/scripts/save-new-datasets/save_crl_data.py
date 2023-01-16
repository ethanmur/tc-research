import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import helper_fns
os.chdir(  "/Users/etmu9498/research/code/scripts/save-new-datasets")
import testing_distance_coords
import find_crl_distance
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts")
import get_p3_heights

def save_all_crl(tcname='all', shift_crl_dist=False, add_dist_coords=False):

    if tcname == 'all':
        tcname_list = [ 'fred', 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tcname]

    for tcname in tcname_list:
        # load the metadata for this tc
        metadata = tc_metadata.all_data( tc= tcname)
        # look at one specific cross section from this tc

        # print( range( len( metadata['dates'])))

        for dataset in range( len( metadata['dates'])):

            # use the function below to process and clean up the tdr data!
            crl_data = save_one_crl( tcname, dataset, metadata, shift_crl_dist, add_dist_coords)

            if crl_data is None:
                return

            # make a name for the dataset
            filename = 'crl-' + metadata[ 'tc_name'].lower() + '-' + metadata[ 'dates'][ dataset] + '-eye-' + metadata[ 'eye_pass'][ dataset] + '.nc'
            # save the newly created crl dataset
            crl_data.to_netcdf('/Users/etmu9498/research/data/crl-new/' + filename)

            print( "New CRL File Created and Saved: " + filename)
    return

# clip the crl dataset to the indices provided in the metadata file, and save it
# this code also preserves the xarray's metadata
# the final xarray dataset is returned through this function
def save_one_crl(tcname, dataset, metadata, shift_crl_dist, add_dist_coords):
    # load the data names for this case
    crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])

    # special 4 ind case: use datasets that have already chopped out microphysical spiral
    # the function that creates these datasets is under save-new-datasets/trim_spiral
    if len( metadata[ 'crl_range'][dataset]) == 4:
        path = "/Users/etmu9498/research/data/crl-spiral-cases" # path to previously created data
    # general 2 ind case: use original dataset
    else:
        path = metadata[ 'crl_path']

    # load the actual crl data for editing
    os.chdir( path)
    crl_data = xr.open_dataset( crl_name)

    # special 4 ind case: use datasets that have already chopped out microphysical spiral
    if len( metadata[ 'crl_range'][dataset]) == 4:
        # i1 and i2 are different in this case because this dataset has already been trimmed down!
        i1 = 0
        i2 = len( crl_data.time)
    # general 2 ind case: use original dataset
    else:
        # define i1 and i2 locally for convenience
        i1 = metadata[ 'crl_range'][dataset][0]
        i2 = metadata[ 'crl_range'][dataset][1]

    # make a copy of the crl data: this will eventually become the new, full dataset!
    crl_new = crl_data.copy()
    # this works really well at dropping all the coords, dims, and vars not listed
    crl_new = crl_new[ ['ProductionDateTime', 'VersionID']]

    # add a new attribute explaining how this dataset has been edited
    disclaimer = 'This Dataset is a subset of the full CRL dataset for ' + metadata['dates'][dataset] + '-2021. Please see the file ' + crl_name + ' for the full dataset.'
    crl_new.attrs[ 'global_att4'] = disclaimer

    disclaimer2 = "Microphysical spirals in the eyes of TC Grace, 8/18, Pass 2 and TC Henri, 8/21, Pass 1 were removed for data continuity purposes. "
    crl_new.attrs[ 'global_att5'] = disclaimer2

    # add time and height as coordinates, not variables!
    cliptime = crl_data.time[ i1:i2]
    crl_new = crl_new.assign_coords( {'time': cliptime })
    crl_new = crl_new.assign_coords( {'H': crl_data.H }) # leaving this name as H to keep consistent with old code


    # find_dist_new_tdr implicitly looks at the current new_crl dataset, so if it's
    # not yet created, it'll blow up :( only run add_dist_coords if the new crl files exist!

    if add_dist_coords:
        # use a helper function to determine distance coords from lat / lon values!
        # new_dist is from tdr based conversions, while new_dist2 is from in situ data
        new_dist = find_crl_distance.find_dist_new_tdr( tcname, dataset)
        print( "TDR Distance axis created")
        new_dist2 = find_crl_distance.find_crl_dist_in_situ( tcname, dataset, returnval='crl')
        print( "In Situ Distance axis created")

        new_dist3 = find_crl_distance.find_crl_dist_psurf( tcname, dataset, returnval='crl')
        print( "Surface Pressure Distance axis created")


        if new_dist2 is None:
            return

        # save the newly created distance array as a coordinate in the crl file!

        crl_new = crl_new.assign_coords( {'tdr_distance': new_dist })
        crl_new.tdr_distance.attrs['long_name'] = 'tdr_distance'
        crl_new.tdr_distance.attrs['units'] = 'meters'
        crl_new.tdr_distance.attrs['description'] = 'Distance from the center of the TC, determined from TDR data'

        # save the in situ distance array in a separate variable
        crl_new = crl_new.assign_coords( {'in_situ_distance': new_dist2 })
        crl_new.in_situ_distance.attrs['long_name'] = 'in_situ_distance'
        crl_new.in_situ_distance.attrs['units'] = 'meters'
        crl_new.in_situ_distance.attrs['description'] = 'Distance from the center of the TC, determined from P-3 height data'

        # saving the surface pressure array
        crl_new = crl_new.assign_coords( {'psurf_distance': new_dist3 })
        crl_new.psurf_distance.attrs['long_name'] = 'psurf_distance'
        crl_new.psurf_distance.attrs['units'] = 'meters'
        crl_new.psurf_distance.attrs['description'] = 'Distance from the center of the TC, determined from in situ surface pressure data'


    # loop through all the variables in crl_data, cropping them to the smaller i1 and i2 values provided
    var_list = list( crl_data.keys() )
    for key_ind in range( len( var_list)):
        # saving 2D arrays
        try:
            # clip the data
            currentdata = crl_data[ var_list[ key_ind]]
            # put everything into numpy form using .data to avoid xr error
            dataclip = currentdata [ i1:i2].data
            # add data to the new, smaller xarray file
            crl_new = crl_new.assign( { var_list[ key_ind]: ( ('time', 'H'), dataclip) }) # used to be ('height':'time')

            # add all the attribute metadata to each variable in the new dataset!
            # get the old metadata
            attr_dict = crl_data[ var_list[ key_ind]].attrs
            # recursively add all keys to new variable!
            for key2 in attr_dict:
                crl_new[ var_list[ key_ind]].attrs[ key2] = attr_dict[ key2]

        # not a 2D array: do nothing
        # all 1D arrays and values have been accounted for here or below
        except (IndexError, ValueError):
            continue

    # add lat and lon manually because it's being a pain to do it automatically
    # add them as coordinates, not variables
    lonclip = crl_data.Lon[ i1:i2].data
    crl_new = crl_new.assign_coords( { 'Lon':  lonclip })
    latclip = crl_data.Lat[ i1:i2].data
    crl_new = crl_new.assign_coords( { 'Lat':  latclip })
    # add metadata for lat and lon
    for key in ['Lon', 'Lat']:
        attr_dict = crl_data[ key].attrs
        # recursively add all keys to new dataset!
        for key2 in attr_dict:
            crl_new[ key].attrs[ key2] = attr_dict[ key2]

    return crl_new









def save_all_new_matrices( tcname='all', add_rmw = False, dist_lim=False):
    if tcname == 'all':
        tcname_list = [ 'fred', 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tcname]

    for tcname in tcname_list:
        # load the metadata for this tc
        metadata = tc_metadata.all_data( tc= tcname)
        # look at one specific cross section from this tc
        for dataset in range( len( metadata['dates'])):
            # use the function below to process and clean up the tdr data!
            crl_data = save_one_p3_height_matrix( tcname, dataset, add_rmw, dist_lim)

            if crl_data is None:
                return

            # make a name for the dataset
            filename = 'crl-' + metadata[ 'tc_name'].lower() + '-' + metadata[ 'dates'][ dataset] + '-eye-' + metadata[ 'eye_pass'][ dataset] + '.nc'
            # save the newly created crl dataset
            crl_data.to_netcdf('/Users/etmu9498/research/data/crl-new-matrices/' + filename)

            print( "New CRL File Created and Saved: " + filename)
    return


def save_one_p3_height_matrix(tcname, dataset, add_rmw, dist_lim=False):
    warnings.filterwarnings("ignore")

    metadata = tc_metadata.all_data( tc= tcname)
    # load the data names for this case
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    # load the actual crl data for editing
    os.chdir( metadata[ 'new_crl_path'])
    crl_data = xr.open_dataset( crl_name)

    # add a new attribute explaining how this dataset has been edited
    disclaimer = ('This Dataset contains both original and height corrected CRL matrices.' +
                ' Results are based on the original file ' + crl_name + '.')
    crl_data.attrs[ 'global_att5'] = disclaimer

    # code to add radius of maximum winds!
    if add_rmw:
        # use a helper function to add rmw data
        new_dist3, new_dist4 = find_crl_distance.find_crl_rmw( tcname, dataset, dist_lim)
        print( "RMW axis created")

        # save the in rmw array in a separate variable
        crl_data = crl_data.assign_coords( {'rmw': new_dist3 })
        crl_data.in_situ_distance.attrs['long_name'] = 'rmw'
        crl_data.in_situ_distance.attrs['units'] = 'unitless'
        crl_data.in_situ_distance.attrs['description'] = 'TC center: RMW = 0. Distance from the TC center to the max in situ wind speed: RMW = 1'

        # the same as rmw above, but with negative values on axes, too!
        crl_data = crl_data.assign_coords( {'rmw_negatives': new_dist4 })
        crl_data.in_situ_distance.attrs['long_name'] = 'rmw_negatives'
        crl_data.in_situ_distance.attrs['units'] = 'unitless'
        crl_data.in_situ_distance.attrs['description'] = 'TC center: RMW = 0. Distance from the TC center to the max in situ wind speed: RMW = 1 or -1'

        '''
        # also save rmw values as floats! this will be convenient for plotting, etc
        crl_data = crl_data.assign_coords( {'rmw': new_dist3 })
        crl_data.in_situ_distance.attrs['long_name'] = 'rmw'
        crl_data.in_situ_distance.attrs['units'] = 'unitless'
        crl_data.in_situ_distance.attrs['description'] = 'TC center: RMW = 0. Distance from the TC center to the max in situ wind speed: RMW = 1'
        '''

    # add new matrices that account for changes in P-3 height!

    # fix matrices before running the helper function to account for height
    # based on code in the make_plots function
    T_2d = crl_data.T.where( crl_data.T.values < 50)

    power_2d = 10 * np.log10( crl_data.P_ch1 )
    power_2d = power_2d.where( power_2d.values > -30)

    wv_2d = crl_data.WVMR.where( crl_data.WVMR.values != 0)
    wv_2d = wv_2d.where( wv_2d.values < 20)

    lsr_2d = crl_data.LSR.where( crl_data.LSR.values < 10)
    lsr_2d = lsr_2d.where( lsr_2d.values > .1)

    # print( 'Original matrices fixed')

    # use a helper fn to get the p-3 height at each position
    p3_heights = get_p3_heights.geth( tcname, dataset)

    # save all p3 height values for each crl timestep
    crl_data = crl_data.assign_coords( {'p3_heights': p3_heights })


    # save the max height to make plotting easier later!
    maxh = np.nanmax( p3_heights) / 1000 # in km

    # getting rid of incorrect ida case
    if maxh > 3.5:
        maxh = 3.5
    crl_data = crl_data.assign_coords( {'H_max': maxh })

    # print("P-3 heights found")

    # use a helper fn to interpolate the data with the new height lims!
    '''
    newh, T_2d = get_p3_heights.interp_data( T_2d, p3_heights)
    newh, power_2d = get_p3_heights.interp_data( power_2d, p3_heights)
    newh, wv_2d = get_p3_heights.interp_data( wv_2d, p3_heights)
    newh, lsr_2d = get_p3_heights.interp_data( lsr_2d, p3_heights)
    '''
    # new function
    # this case accounts for the annoying extra 1 km of noise in grace 8/18 data ;/
    # fix this line again if adding 8/16 data too!
    if tcname == 'grace' and dataset < 3:
        newh, T_2d = get_p3_heights.interp_data2( T_2d, crl_data.H, p3_heights, grace_case=17)
        newh, power_2d = get_p3_heights.interp_data2( power_2d, crl_data.H, p3_heights, grace_case=17)
        newh, wv_2d = get_p3_heights.interp_data2( wv_2d, crl_data.H, p3_heights, grace_case=17)
        newh, lsr_2d = get_p3_heights.interp_data2( lsr_2d, crl_data.H, p3_heights, grace_case=17)

    elif tcname == 'grace' and dataset > 2 and dataset < 6:
        newh, T_2d = get_p3_heights.interp_data2( T_2d, crl_data.H, p3_heights, grace_case=18)
        newh, power_2d = get_p3_heights.interp_data2( power_2d, crl_data.H, p3_heights, grace_case=18)
        newh, wv_2d = get_p3_heights.interp_data2( wv_2d, crl_data.H, p3_heights, grace_case=18)
        newh, lsr_2d = get_p3_heights.interp_data2( lsr_2d, crl_data.H, p3_heights, grace_case=18)

    else:
        newh, T_2d = get_p3_heights.interp_data2( T_2d, crl_data.H, p3_heights)
        newh, power_2d = get_p3_heights.interp_data2( power_2d, crl_data.H, p3_heights)
        newh, wv_2d = get_p3_heights.interp_data2( wv_2d, crl_data.H, p3_heights)
        newh, lsr_2d = get_p3_heights.interp_data2( lsr_2d, crl_data.H, p3_heights)

    # print( "Data Interpolated")

    # calculate relative humidity and temp anomaly with corrected heights and matrices...
    # there would be errors without these corrections!

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


    # make a new set of vertical coordinates for new matrices
    crl_data = crl_data.assign_coords( {'H_new': newh })

    # add the new matrices to the xarray dataset!
    # does flipping time and height work? or break things?
    # print( type( T_2d))
    # print( type( crl_rh))
    # print( np.shape( T_2d))
    # print( np.shape( crl_rh) )

    '''
    print( type( xr.DataArray( lsr_2d)))
    print( np.shape( xr.DataArray( lsr_2d)))
    print( type( crl_rh))
    print( np.shape( crl_rh))
    '''

    # crl_data = crl_data.assign( { 'rh_new': ( ('time', 'H_new'), crl_rh.data) })

    crl_data = crl_data.assign( { 'T_new': xr.DataArray( T_2d)}) # ( ('height', 'time'), T_2d.transpose() ) })
    crl_data = crl_data.assign( { 'power_new': xr.DataArray( power_2d)}) # ( ('H_new', 'time'), power_2d) })

    crl_data = crl_data.assign( { 'rh_new': crl_rh })
    crl_data = crl_data.assign( { 'T_anomaly': xr.DataArray( temp_anomaly_obs )})

    crl_data = crl_data.assign( { 'wvmr_new': xr.DataArray( wv_2d)}) # ( ('H_new', 'time'), wv_2d) })
    crl_data = crl_data.assign( { 'lsr_new': xr.DataArray( lsr_2d)}) # ( ('H_new', 'time'), lsr_2d) })

    '''
    print( type( xr.DataArray( lsr_2d)))
    print( np.shape( xr.DataArray( lsr_2d)))
    print( type( crl_rh.values))
    print( np.shape( crl_rh))
    '''

    warnings.filterwarnings("default")

    return crl_data
