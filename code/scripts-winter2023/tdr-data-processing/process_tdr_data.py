import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
os.chdir("/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023

def save_all_tdr(tcname='all'):
    crl_data_root = "/Users/etmu9498/research/data/tdr-original/"
    yearlist, namelist, filelist = helper_fns_winter2023.get_tdr_datasets( tcname)

    for yeari, yearval  in enumerate( yearlist):
        for tci, tcval in enumerate( namelist[yeari]):
            for filei, fileval in enumerate( filelist[yeari][tci][0]):
                # print("year " + yearval + ", tc " + tcval + ", file " + fileval)

                # inbound case: actually combine the datasets!
                # skip outbound cases, they've already been added during inbound cases
                if 'inbound' in fileval:
                    inbound_name = fileval
                    outbound_name = fileval[:-11] + "_outbound.nc"
                    # make sure that the respective outbound file exists, too
                    if outbound_name in filelist[yeari][tci][0]:

                        # use the function below to process and clean up the tdr data!
                        os.chdir("/Users/etmu9498/research/data/tdr-original/"+yearval+'/'+tcval + '/')
                        tdr_data = save_one_tdr( tcname, inbound_name, outbound_name)

                        # make a name for the dataset
                        date = inbound_name[0:13]
                        filename = 'tdr_combined_' +tcval+ '_' +date+ '.nc'
                        
                        # make correct folder dependencies
                        # see if there's already a tc year folder availible
                        os.chdir("/Users/etmu9498/research/data/tdr-one-pass/")
                        output_folder = yearval
                        if not os.path.isdir( output_folder):
                            os.makedirs( output_folder)
                            print( 'New folder created: ' + output_folder)
                        # see if there's already a tc name folder availible
                        os.chdir("/Users/etmu9498/research/data/tdr-one-pass/" + yearval )
                        output_folder = tcval
                        if not os.path.isdir( output_folder):
                            os.makedirs( output_folder)
                            print( 'New folder created: ' + output_folder)

                        # save the newly created tdr dataset
                        tdr_data.to_netcdf('/Users/etmu9498/research/data/tdr-one-pass/' + yearval + "/" + tcval + "/" + filename)
                    else:
                        print("Error: " + outbound_name + " not found in folder of interest. No combined TDR file will be created.")

            print( "TDR Files Created and Saved for TC " + tcval)
    return
    

# combine one inbound tdr dataset with the correct outbound dataset
# this code also preserves the xarray's metadata
# the final xarray dataset is returned through this function for saving in save_all_tdr()
def save_one_tdr(tcname, inbound_name, outbound_name):
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # make a copy of the inbound data: this will eventually become the new, full dataset!
    tdr_data = inbound_data.copy()
    # this works really well at dropping all the coords, dims, and vars not listed!!
    tdr_data = tdr_data[ ['time', 'height', 'heading']]

    # create a new distance array, joining the two radius arrays
    # convert from xr to np to concatonate arrays
    # which dataset should get the - sign??? I think inbound works better
    x_in = - inbound_data.radius
    x_out = outbound_data.radius
    x_in = xr.DataArray.to_numpy( x_in)
    x_out = xr.DataArray.to_numpy( x_out)
    # flip outbound data so everything is in order when concatonating!
    x_out = np.flip( x_out, 0)
    # create one distance array! add a 0 km data point to keep things even, too :)
    x_in = np.append( [0.0], x_in)
    distance_array = np.append( x_out, x_in)

    # add the distance array as a coordinate to the new dataset
    tdr_data = tdr_data.assign_coords( {'distance': distance_array})
    # add metadata to the distance coord
    tdr_data.distance.attrs['UNITS'] = 'meters'
    tdr_data.distance.attrs['STANDARD_NAME'] = 'distance'
    tdr_data.distance.attrs['AXIS'] = 'X'

    # last step: loop through the tdr dataset to add parameters, all in the proper format!

    # save the names of all the variables in the original xarray dataset to:
    var_list = list( inbound_data.keys() )

    # i got rid of the code (below) to add new names because it was messing with compatability with my older code
    # make a list of new, clearer names to save the variables to!
    new_names = ['azimuth', 'altitude', 'latitude', 'longitude', 'seconds', 'u_air', 'v_air', 'w_air',
                 'vgw_air', 'ws_air', 'radial_wind', 'tangential_wind', 'vertical_wind', 'wind_speed', 'reflectivity']

    # save all the new values to the array!
    for key_ind in range( len( var_list)):
        # saving 2D arrays
        try:
            # select the correct data, and flip the outbound data to match the flipped outbound distance!
            out_vals = outbound_data[ var_list[ key_ind]].isel(time=0).isel(heading=0).transpose().values
            out_vals = np.flip( out_vals, 1)
            in_vals = inbound_data[ var_list[ key_ind]].isel(time=0).isel(heading=0).transpose().values

            # special case for velocities: need to flip sign of inbound data for correct plotting later!!!
            if var_list[ key_ind] == 'Radial_wind' or var_list[ key_ind] == "Vertical_wind":
                in_vals = in_vals # * -1
                # print( "vel case!")

            # add one empty vertical column of data to account for empty 0 km data point
            nan_pad = np.empty(  (len( outbound_data.height), 1) ) # create an empty array
            nan_pad[:] = np.nan # fill it with nans

            # combine the three new datasetS!
            in_vals = np.append( nan_pad, in_vals, axis=1)
            all_vals = np.append( out_vals, in_vals, axis=1)

            tdr_data = tdr_data.assign( { var_list[ key_ind]: ( ('height', 'distance'), all_vals) })
            # code to add new names to tdr variables
            # tdr_data = tdr_data.assign( { new_names[ key_ind]: ( ('height', 'distance'), all_vals) })
        # saving 1D arrays
        except (IndexError, ValueError):
            # set up the correct variables
            in_vals = inbound_data[ var_list[ key_ind]]
            out_vals = outbound_data[ var_list[ key_ind]]
            in_vals = xr.DataArray.to_numpy( in_vals)
            out_vals = xr.DataArray.to_numpy( out_vals)

            # flip outbound data so everything is in order when concatonating
            out_vals = np.flip( out_vals, 0)
            in_vals = np.append( [np.nan], in_vals) # use nan, not 0, so no data shows up here
            all_vals = np.append( out_vals, in_vals)

            tdr_data = tdr_data.assign( { var_list[ key_ind]: (('distance'), all_vals ) })
            # tdr_data = tdr_data.assign( { new_names[ key_ind]: (('distance'), all_vals ) })

        # add all the attribute metadata to the new dataset!
        # get the old metadata
        attr_dict = inbound_data[ var_list[ key_ind]].attrs
        # recursively add all keys to new dataset
        for key2 in attr_dict:
            tdr_data[ var_list[ key_ind]].attrs[ key2] = attr_dict[ key2]
            # code to add new names to variables
            # tdr_data[ new_names[ key_ind]].attrs[ key2] = attr_dict[ key2]

    return tdr_data
