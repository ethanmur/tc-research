# This saving process is very similar to the save_crl_data.py script found in the same
# folder!
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
import sys
sys.path.append("/Users/etmu9498/research/code/scripts/in-situ-scripts/")
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts/")
import load_in_situ_data

def save_all_in_situ(tcname='all', shift_crl_dist=False, add_dist_coords=False):

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
            in_situ_data = save_one_in_situ( tcname, dataset, metadata, shift_crl_dist, add_dist_coords)

            if in_situ_data is None:
                return

            # make a name for the dataset
            filename = 'in-situ-' + metadata[ 'tc_name'].lower() + '-' + metadata[ 'dates'][ dataset] + '-eye-' + metadata[ 'eye_pass'][ dataset] + '.nc'
            # save the newly created crl dataset
            in_situ_data.to_netcdf('/Users/etmu9498/research/data/in-situ-new/' + filename)

            print( "New In Situ File Created and Saved: " + filename)
    return

def save_one_in_situ(tcname, dataset, metadata, shift_crl_dist, add_dist_coords):
    # load the data names for this case
    insitu_list = make_plots.load_flight_level( metadata['in_situ_path'], print_files=False)
    insitu_name = tc_metadata.choose_in_situ_date( metadata[ 'dates'][ dataset], insitu_list)

    # load the actual data for editing
    in_situ_data = load_in_situ_data.load_in_situ( metadata['in_situ_path'], insitu_name)
    newData = in_situ_data.copy()

    print( 'converted to nc file')

    # use a helper script to find the new distance array for the in situ data!
    distance = find_crl_distance.find_crl_dist_in_situ( tcname, dataset, returnval='in-situ')

    print( 'in situ distance found')

    # save the in situ distance array in a separate variable
    # newData= newData.assign_coords( {'distance': distance })
    # newData.distance.attrs['long_name'] = 'distance'
    # newData.distance.attrs['units'] = 'meters'
    # newData.distance.attrs['description'] = 'Distance from the center of the TC, determined from minimum P-3 height during an eye pass'

    # also save the distance as a variable for convenience!
    newData = newData.assign( {'distance': (('time'), distance ) })
    newData.distance.attrs['long_name'] = 'distance'
    newData.distance.attrs['units'] = 'meters'
    newData.distance.attrs['description'] = 'Distance from the center of the TC, determined from minimum P-3 height during an eye pass'

    return newData
