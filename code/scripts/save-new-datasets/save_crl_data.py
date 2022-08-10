import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import helper_fns
os.chdir(  "/Users/etmu9498/research/code/scripts/save-new-datasets")
import testing_distance_coords
import find_crl_distance

def save_all_crl(tcname='all', shift_crl_dist=False, add_dist_coords=False):

    if tcname == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tcname]

    for tcname in tcname_list:
        # load the metadata for this tc
        metadata = tc_metadata.all_data( tc= tcname)
        # look at one specific cross section from this tc
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
    # load the actual crl data for editing
    os.chdir( metadata[ 'crl_path'])
    crl_data = xr.open_dataset( crl_name)

    # make a copy of the crl data: this will eventually become the new, full dataset!
    crl_new = crl_data.copy()
    # this works really well at dropping all the coords, dims, and vars not listed!!
    crl_new = crl_new[ ['ProductionDateTime', 'VersionID']]

    # add a new attribute explaining how this dataset has been edited
    disclaimer = 'This Dataset is a subset of the full CRL dataset for ' + metadata['dates'][dataset] + '-2021. Please see the file ' + crl_name + ' for the full dataset.'
    crl_new.attrs[ 'global_att4'] = disclaimer

    # define i1 and i2 locally for convenience
    i1 = metadata[ 'crl_range'][dataset][0]

    # special 4 ind case to chop out
    if len( metadata[ 'crl_range'][dataset]) == 4:
        i2 = metadata[ 'crl_range'][dataset][3]
    # general 2 ind case
    else:
        i2 = metadata[ 'crl_range'][dataset][1]

    # add time and height as coordinates, not variables!
    cliptime = crl_data.time[ i1 : i2 ]
    crl_new = crl_new.assign_coords( {'time': cliptime })
    crl_new = crl_new.assign_coords( {'H': crl_data.H }) # leaving this name as H to keep consistent with old code

    # you might have to comment out the next four lines of code to get things to work!
    # find_dist_new_tdr implicitly looks at the current new_crl dataset, so if it's
    # not yet created, it'll blow up :( only run add_dist_coords if the files exist!

    if add_dist_coords:
        # use a helper function to determine distance coords from lat / lon values!
        # new_dist is from tdr based conversions, while new_dist2 is from in situ data
        new_dist = find_crl_distance.find_dist_new_tdr( tcname, dataset)
        new_dist2 = find_crl_distance.find_crl_dist_in_situ( tcname, dataset, returnval='crl')

        if new_dist2 is None:
            return

        # shifting the crl center point distances using metadata stored in all_plots
        # this is optional and defaults to False
        if shift_crl_dist:
            shift = metadata[ 'shift'][ dataset]
            stretch = metadata[ 'stretch'][ dataset]
            new_dist = new_dist * stretch + shift

        # save the newly created distance array as a coordinate in the crl file!
        crl_new = crl_new.assign_coords( {'tdr_distance': new_dist })
        crl_new.tdr_distance.attrs['long_name'] = 'tdr_distance'
        crl_new.tdr_distance.attrs['units'] = 'meters'
        crl_new.tdr_distance.attrs['description'] = 'Distance from the center of the TC, determined from TDR data'

        # save the in situ distance array in a separate variable
        crl_new = crl_new.assign_coords( {'in_situ_distance': new_dist2 })
        crl_new.in_situ_distance.attrs['long_name'] = 'in_situ_distance'
        crl_new.in_situ_distance.attrs['units'] = 'meters'
        crl_new.in_situ_distance.attrs['description'] = 'Distance from the center of the TC, determined from in situ P-3 data'

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
            crl_new = crl_new.assign( { var_list[ key_ind]: ( ('height', 'time'), dataclip) })

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
