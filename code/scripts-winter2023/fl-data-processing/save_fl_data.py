# This saving process is very similar to the save_crl_data.py script found in the same
# folder!
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
import sys
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/fl-data-processing")
import find_fl_dists_rmws


# this function is the main loop to locally save multiple new flight level datasets.
# inputs:
############
# tc: can either be 'all' to do this for all years, a year string like
#     '2021' or '2022', or a dictionary with years as the key and files to look at in
#     a list as the values. The last option is the most useful for smaller tests!
# add_dist_coords: a dictionary of possible distance coordinates to add / other optional
#     inputs. All the inputs default to False (kinda useless lol but prevents errors!)
#     'vars_to_save' is a bit different. If = 'default', the default variables included in
#     the save_one_fl() function are used. otherwise, you can pass a list of variables through
#     this field for custom behavior!
def save_tcs( tc='all', add_dist_coords={'center_dist': False, 'psurf_dist': False, 'rmw': False, 'vars_to_save': 'default'}):

    data_path = "/Users/etmu9498/research/data/"
    fl_data_path = "/Users/etmu9498/research/data/in-situ-noaa-full/"

    # do this for all flight level datasets
    if tc == 'all':
        # make a list of years where crl data is present
        os.chdir( data_path)
        yearlist = [name for name in os.listdir('in-situ-noaa-full')
            if os.path.isdir(os.path.join('in-situ-noaa-full', name))]

        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        for yeari, yearval in enumerate( yearlist):
            filelist.append( make_plots.load_flight_level( fl_data_path + yearval, print_files=False) )

    # do this for just one year
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):
        yearlist = [ tc]
        filelist = [ make_plots.load_flight_level( fl_data_path + tc, print_files=False)]

    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")

    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be processed: " + str( filecount))


    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            # use the function below to process the crl data separately!
            fl_data = save_one_fl( yearval, fileval, add_dist_coords= add_dist_coords)

            # make a name for the dataset
            filename = filelist[ yeari][ filei] [ : -3] + "_processed.nc"
            # save the newly created crl dataset

            # make sure the folder exists!
            os.chdir("/Users/etmu9498/research/data/in-situ-noaa-processed")
            output_folder = yearval
            if not os.path.isdir( output_folder):
                os.makedirs( output_folder)
                print( 'New folder created: ' + output_folder)
            #else:
            #    print( 'Folder ' + output_folder + ' already exists')

            # go to the new folder and save the data!
            fl_data.to_netcdf('/Users/etmu9498/research/data/in-situ-noaa-processed/' + yearval + "/" + filename)
            print( "New In Situ File Created and Saved: " + filename + "\n")
    return







# trim down and create new time, distance, and rmw axes for one flight level dataset!
# inputs:
#########
# yearval: the year of this tc occurence
# fileval: the name of this tc dataset (in the format of "20210929H2_sam.nc")
# return:
#########
# the edited fl dataset.

def save_one_fl( yearval, fileval, add_dist_coords={'center_dist': False, 'psurf_dist': False, 'rmw': False, 'vars_to_save': 'default'}):

    # use the default coordinates case
    if add_dist_coords['vars_to_save'] == 'default':
        # 'float_time', 'Time',
        input_vars = [ 'HT.d', 'THETA.d', 'THETAV.d', 'WS.d', 'UWZ.d',
                'SfmrRainRate.1', 'THETAE.d', 'MR.d', 'TA.d', 'PSURF.d', 'LATref', 'LONref'] # change me!

    # use an input list of input_vars! This takes some trust for the user; make sure all the vars are applicable!
    elif type( add_dist_coords['vars_to_save']) == type( []):
        input_vars = add_dist_coords['vars_to_save']

    # load the raw, original data!
    fl_new_path = "/Users/etmu9498/research/data/in-situ-noaa-full/" + yearval
    os.chdir(fl_new_path)
    fl_data = xr.open_dataset( fileval, decode_times=False)
    # since the tc names are nicely saved at the same spot on all datasets, we can pull
    # just their names pretty easily!
    tcname = fileval[ 11 : -3]


    # make a copy of the crl data: this will eventually become the new, full dataset!
    # fl_new = fl_data.copy()


    # creating the time (decimal) axis
    # interval_str holds the start and end times as a string. cut down interval_str to get the start hour, min, sec!
    interval_str = fl_data.attrs['TimeInterval']
    h, m, s = float( interval_str[0:2]), float( interval_str[3:5]), float( interval_str[6:8])
    start_time = h + m / 60 + s / 3600

    # create the time array manually. This is possible because time values increase consistently at 1 second intervals!
    time = np.empty( ( len( fl_data['Time'])))
    for timei in range( len( fl_data['Time'])):
        # add to time array
        time[ timei] = start_time + timei / 3600

    # fl_new['float_time'] = time
    # trim down the dataset to just the axes of interest!
    # fl_new = fl_new[ input_vars]
    dvs = {}
    units = []
    description = []
    sample = []
    output = []
    valid = []


    # optional code: find the radial distance from the tc center!
    # this takes ~ 30 seconds to run per case, and a valid flight track must be
    # provided for this to work!
    if add_dist_coords['center_dist'] or add_dist_coords['rmw']:

        # use a helper function to find the distance axis
        dists, rmw = find_fl_dists_rmws.find_center_dist_rmw( yearval, fileval, tcname, fl_data, time)

        # final step: add the distance and or rmw axes to the dataset!
        if add_dist_coords['rmw']:
            # fl_new['rmw'] = rmw
            dvs.update( { 'rmw': ('time', rmw) })
            units.append( 'Unitless' )
            description.append( 'Normalized distance to the radius of maximum winds.' )
            sample.append( 1.0 )
            output.append( 1.0 )
            valid.append( [0.0, 200.0] )

        if add_dist_coords['center_dist']:
            # fl_new['center_dist'] = dists
            dvs.update( { 'center_dist': ('time', dists) })
            units.append( 'Km' )
            description.append( 'Distance from the TC center.' )
            sample.append( 1.0 )
            output.append( 1.0 )
            valid.append( [0.0, 2000.0] )

    # add new fields to the dataset!
    dvs = {}
    for fieldi, fieldval in enumerate( input_vars):
        field = fl_data[ fieldval].values
        # add this field to the dataframe!
        dvs.update( { fieldval: ('time', field) })


    # this code hasn't been implemented yet
    if add_dist_coords['psurf_dist']:
        pass

    # next step: create the dataset and update units!
    fl_new = xr.Dataset( data_vars=dvs, coords={'time': time})
    for fieldi, fieldval in enumerate( fl_new.variables):
        if fieldval in fl_data.variables:
            attribs = fl_data[ fieldval].attrs
            # add units, etc metadata present in the original flight level dataset!
            if 'units' in attribs:
                fl_new[fieldval].attrs['units'] = fl_data[ fieldval].units
            if 'Description' in attribs:
                # description.append( fl_data[ fieldval].Description )
                fl_new[fieldval].attrs['Description'] = fl_data[ fieldval].Description
            if 'SampleRate' in attribs:
                # sample.append( fl_data[ fieldval].SampleRate )
                fl_new[fieldval].attrs['SampleRate'] = fl_data[ fieldval].SampleRate
            if 'OutputRate' in attribs:
                # output.append( fl_data[ fieldval].OutputRate )
                fl_new[fieldval].attrs['OutputRate'] = fl_data[ fieldval].OutputRate
            if 'ValidRange' in attribs:
                # valid.append( fl_data[ fieldval].ValidRange )
                fl_new[fieldval].attrs['ValidRange'] = fl_data[ fieldval].ValidRange
            if 'Min' in attribs:
                fl_new[fieldval].attrs['Min'] = fl_data[ fieldval].Min
            if 'Max' in attribs:
                fl_new[fieldval].attrs['Max'] = fl_data[ fieldval].Max

    # add a new attribute explaining how this dataset has been edited
    author = 'Created by NOAA HRD. Link to data: https://www.aoml.noaa.gov/ftp/hrd/data/flightlevel/'
    fl_new.attrs['Author'] = author
    disclaimer = 'Edited by: Ethan Murray (etmu9498@colorado.edu)'
    fl_new.attrs[ 'Editor'] = disclaimer
    disclaimer = 'This Dataset is based off the dataset ' + fileval + ". Please see that file for the original data."
    fl_new.attrs[ 'Attribution'] = disclaimer


    return fl_new
