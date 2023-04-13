# Author: Ethan Murray
# Created:
# Edited: 4/3/23
# Automatically trim and save all relevant flight level datasets.
# The saving process is very similar to the save_crl_data.py script found in the scripts-winter2023 folder

# import...
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import xarray as xr
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/fl-data-processing")
import find_fl_dists_rmws
import correct_fl_data_fields


# this function is the main loop to locally save multiple new flight level datasets.
# inputs:
############
# tc: can either be 'all' to do this for all years, a year string like
#     '2021' or '2022', or a dictionary with years as the key and files to look at in
#     a list as the values. The last option is the most useful for smaller tests!
# add_dist_coords: a dictionary of possible distance coordinates to add / other optional
#     inputs. All the inputs default to False.
#     'vars_to_save' is a bit different. If = 'default', the default variables included in
#     the save_one_fl() function are used. otherwise, you can pass a list of variables through
#     this field for different behavior!
# correct_data tells the script whether or not to process the flight level data
#     (removing too high / low legs, anomalously high wv, etc).
#     You can give it either 'full' to not process the data at all (keep the full dataset),
#     'correct' to just correct the issue data points, or 'trim' to do the corrections and delete problem data. 
#     'correct' also leaves the relative vorticity as it is: it doesn't apply any filters to it!
def save_tcs( tc='all', add_dist_coords={'center_dist': False, 'psurf_dist': False, 'rmw': False, 'vars_to_save': 'default'}, correct_data='correct'):

    data_path = "/Users/etmu9498/research/data/"
    fl_data_path = "/Users/etmu9498/research/data/in-situ-noaa-full/"

    # case 1: do this for all flight level datasets
    if tc == 'all':
        # make a list of years where crl data is present
        os.chdir( data_path)
        yearlist = [name for name in os.listdir('in-situ-noaa-full')
            if os.path.isdir(os.path.join('in-situ-noaa-full', name))]
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        for yeari, yearval in enumerate( yearlist):
            filelist.append( helper_fns.load_flight_level( fl_data_path + yearval, print_files=False) )

    # case 2: do this for just one year
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):
        yearlist = [ tc]
        filelist = [ helper_fns.load_flight_level( fl_data_path + tc, print_files=False)]

    # case 3: do this for a specific dictionary of files:
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
            fl_data = save_one_fl( yearval, fileval, add_dist_coords= add_dist_coords, correct_data=correct_data)

            # final goal: save the newly created crl dataset
            # make a name for the dataset
            filename = filelist[ yeari][ filei] [ : -3] + "_processed.nc"

            if correct_data == 'correct':
                # make sure the folder exists!
                path = "/Users/etmu9498/research/data/in-situ-noaa-processed/"
            elif correct_data == 'full':
                path = "/Users/etmu9498/research/data/in-situ-noaa-not-processed/"
            elif correct_data == 'trim':
                path = "/Users/etmu9498/research/data/in-situ-noaa-trimmed/"

            os.chdir(path)
            output_folder = yearval
            if not os.path.isdir( output_folder):
                os.makedirs( output_folder)
                print( 'New folder created: ' + output_folder)
            # go to the new folder and save the data!
            fl_data.to_netcdf(path + yearval + "/" + filename)
            print( "New In Situ File Created and Saved: " + filename + "\n")
    return



# Create one processed flight level dataset by:
# trimming down the saved variables to a specific list,
# and creating new time, distance, and rmw axes, if applicable.
# inputs:
#########
# yearval: the year of this tc occurence
# fileval: the name of this tc dataset (in the format of "20210929H2_sam.nc")
# return: a new trimmed netcdf file with the additional coordinates.
#########
def save_one_fl( yearval, fileval, add_dist_coords={'center_dist': False, 'psurf_dist': False, 'rmw': False, 'vars_to_save': 'default'}, correct_data='correct', limits=False, inputdata=False):
    # option 1: use the default variables case
    if add_dist_coords['vars_to_save'] == 'default':
        # 'float_time', 'Time',
        input_vars = [ 'HT.d', 'THETA.d', 'THETAV.d', 'WS.d', 'WD.d', 'UWZ.d', 'UWX.d', 'UWY.d', 'HUM_REL.d',
                'SfmrRainRate.1', 'THETAE.d', 'MR.d', 'TA.d', 'PSURF.d', 'LATref', 'LONref'] # change me!

    # option 2: use an input list of input_vars! This takes some trust for the user; manually check that all the vars are applicable!
    elif type( add_dist_coords['vars_to_save']) == type( []):
        input_vars = add_dist_coords['vars_to_save']

    # special case: use provided (often previously processed) flight level data!
    if inputdata:
        fl_data = inputdata
    # default case:
    # load the raw, original flight level data!
    else:
        fl_new_path = "/Users/etmu9498/research/data/in-situ-noaa-full/" + yearval
        os.chdir(fl_new_path)
        fl_data = xr.open_dataset( fileval, decode_times=False)
        tcname = fileval[ 11 : -3]

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

    # add new fields to the dataset!
    dvs = {}

    # save the original p3 height here... it was being iterated over and edited in an
    # unhelpful way!
    p3_height_new = np.copy( fl_data[ 'HT.d'])

    # correct the p-3 height array outside of the loop... will impact all other measurements!
    bad_inds = np.union1d( np.where( fl_data[ 'HT.d'].values < 2500.)[0], np.where( fl_data[ 'HT.d'].values > 4000.)[0] )
    # the correct inds are all inds - bad_inds!
    correct_inds = np.setdiff1d( np.arange( 0, len( fl_data[ 'HT.d'].values)), bad_inds)
    # fl_data[ 'HT.d'][ bad_inds] = np.nan
    p3_height_new[ bad_inds] = np.nan

    print(np.shape(p3_height_new))
    print(len(p3_height_new[np.where(~ np.isnan( p3_height_new))]))
    print(len(fl_data['HT.d'][np.where(~ np.isnan(fl_data['HT.d']))]))
    print(type(p3_height_new))
    print(p3_height_new)

    # make a field list for adding things in a separate loop... doing it in the same loop was causing problems :/
    field_list = []


    #############
    ## 3/17 new code: account for special cases of creating dataset when limits are applied!
    #############
    if limits:
        # fieldval is the NAME of the dataset, not the actual data...
        for fieldi, fieldval in enumerate( input_vars):
            
            # better case: correct for inconsistencies in the flight level data! will improve stats later
            # DON'T correct the height field... leads to issues later!
            if correct_data == 'trim': # and fieldval != 'HT.d':
                # give the helper script the data field of interest, the p-3 height array (too high or low is an issue),
                # and the time array (needed?? maybe not)
                field = correct_fl_data_fields.one_field_trim( fieldval, fl_data[ fieldval].values, p3_height_new, bad_inds, time, fileval, fl_data[ 'TA.d'].values)
            
            elif correct_data == 'correct':
                field = correct_fl_data_fields.one_field_correct( fieldval, fl_data[ fieldval].values, p3_height_new, time, fileval, fl_data[ 'TA.d'].values)

            # original case: just add the raw fl data fields to the array!
            elif correct_data == 'full':
                field = fl_data[ fieldval].values
            field_list.append( field[ limits[0] : limits[1] ])
    
    # original, normal case:
    else:
        # fieldval is the NAME of the dataset, not the actual data...
        for fieldi, fieldval in enumerate( input_vars):

            # better case: correct for inconsistencies in the flight level data! will improve stats later
            # DON'T correct the height field... leads to issues later!
            if correct_data == 'trim': # and fieldval != 'HT.d':
                # give the helper script the data field of interest, the p-3 height array (too high or low is an issue),
                # and the time array (needed?? maybe not)
                field = correct_fl_data_fields.one_field_trim( fieldval, fl_data[ fieldval].values, p3_height_new, bad_inds, time, fileval, fl_data[ 'TA.d'].values)
            
            elif correct_data == 'correct':
                field = correct_fl_data_fields.one_field_correct( fieldval, fl_data[ fieldval].values, p3_height_new, time, fileval, fl_data[ 'TA.d'].values)

            # original case: just add the raw fl data fields to the array!
            elif correct_data == 'full':
                field = fl_data[ fieldval].values

            field_list.append( field)


    for fieldi, fieldval in enumerate( field_list):
        # add this field to the dataframe!
        dvs.update( { input_vars[ fieldi]: ('time', fieldval) })

    # optional code: find the radial distance from the tc center!
    # this takes ~ 30 seconds to run per case, and a valid flight track must be
    # provided for this to work.
    if add_dist_coords['center_dist'] or add_dist_coords['rmw']:
        # use a helper function to find the distance axis
        dists, rmw = find_fl_dists_rmws.find_center_dist_rmw( yearval, fileval, tcname, fl_data, time)

        # add these fields to the original dataset for easy access below!
        fl_data['center_dist'] = dists
        fl_data['rmw'] = rmw

        # final step: add the distance and or rmw axes to the dataset!
        if add_dist_coords['rmw']:
            dvs.update( { 'rmw': ('time', rmw) })
        if add_dist_coords['center_dist']:
            dvs.update( { 'center_dist': ('time', dists) })

    # this code hasn't been implemented yet
    if add_dist_coords['psurf_dist']:
        pass


    # add relative vorticity to the dataset!
    # using the method found in kossin and eastin 2001
    # give the helper fn total wind speed in m/s and radial distance in meters
    # return processed / smoothed relative vorticity in units of 10^-4 s-1

    # spectial 'full' data case: correct vorticity before adding it to the array
    if correct_data=='full':
        if 'WS.d' in fl_data.keys() and 'center_dist' in fl_data.keys():
            wd = 60 # 14
            avg = False
            rel_vort = find_fl_dists_rmws.calc_rel_vort( fl_data['WS.d'].values, fl_data['center_dist'].values * 1000, double_average=avg, window=wd)
            # add it to the dataframe for nicer loops!
            dvs.update( { 'rel_vort': ('time', rel_vort) })
        else:
            print( "relative vorticity not added")
    # don't do any smoothing!
    else:
        if 'WS.d' in fl_data.keys() and 'center_dist' in fl_data.keys():
            wd = 1
            avg = False
            rel_vort = find_fl_dists_rmws.calc_rel_vort( fl_data['WS.d'].values, fl_data['center_dist'].values * 1000, double_average=avg, window=wd)
            # add it to the dataframe for nicer loops!
            dvs.update( { 'rel_vort': ('time', rel_vort) })
        else:
            print( "relative vorticity not added")


    # next step: create the dataset and update units for existing variables!
    fl_new = xr.Dataset( data_vars=dvs, coords={'time': time})


    for fieldi, fieldval in enumerate( fl_new.variables):
        if fieldval in fl_data.variables:
            attribs = fl_data[ fieldval].attrs
            # add units, etc metadata present in the original flight level dataset!
            if 'units' in attribs:
                fl_new[fieldval].attrs['units'] = fl_data[ fieldval].units
            if 'Description' in attribs:
                fl_new[fieldval].attrs['Description'] = fl_data[ fieldval].Description
            if 'SampleRate' in attribs:
                fl_new[fieldval].attrs['SampleRate'] = fl_data[ fieldval].SampleRate
            if 'OutputRate' in attribs:
                fl_new[fieldval].attrs['OutputRate'] = fl_data[ fieldval].OutputRate
            if 'ValidRange' in attribs:
                fl_new[fieldval].attrs['ValidRange'] = fl_data[ fieldval].ValidRange
            if 'Min' in attribs:
                fl_new[fieldval].attrs['Min'] = fl_data[ fieldval].Min
            if 'Max' in attribs:
                fl_new[fieldval].attrs['Max'] = fl_data[ fieldval].Max

    # manually add units for rmw and center_dist variables!
    fl_new['center_dist'].attrs['units'] = 'Km'
    fl_new['center_dist'].attrs['Description'] = 'Distance from the TC center.'
    fl_new['center_dist'].attrs['SampleRate'] = 1.0
    fl_new['center_dist'].attrs['OutputRate'] = 1.0
    fl_new['center_dist'].attrs['ValidRange'] = [0.0, 2000.0]

    fl_new['rmw'].attrs['Units'] = 'Unitless'
    fl_new['rmw'].attrs['Description'] = 'Normalized distance to the radius of maximum winds'
    fl_new['rmw'].attrs['SampleRate'] = 1.0
    fl_new['rmw'].attrs['OutputRate'] = 1.0
    fl_new['rmw'].attrs['ValidRange'] = [0.0, 200.0]


    # add a new attribute explaining how this dataset has been edited
    author = 'Created by NOAA HRD. Link to data: https://www.aoml.noaa.gov/ftp/hrd/data/flightlevel/'
    fl_new.attrs['Author'] = author
    disclaimer = 'Edited by: Ethan Murray (etmu9498@colorado.edu)'
    fl_new.attrs[ 'Editor'] = disclaimer
    disclaimer = 'This Dataset is based off the dataset ' + fileval + ". Please see that file for the original data."
    fl_new.attrs[ 'Attribution'] = disclaimer

    if 'WS.d' in fl_data.keys() and 'center_dist' in fl_data.keys():
        # adding units for relative vorticity
        fl_new['rel_vort'].attrs['Units'] = '* 10^-4 s^-1'
        fl_new['rel_vort'].attrs['Description'] = '1D relative vorticity. Calculated using the equation found on page 1083 of Kossin and Eastin 2001.'
        if avg:
            avgtext = "two times."
        else:
            avgtext = 'one time.'
        fl_new['rel_vort'].attrs['Description2'] = 'Window for averaging = ' + str(wd) + '. Averaged ' + avgtext
        fl_new['rel_vort'].attrs['SampleRate'] = 1.0
        fl_new['rel_vort'].attrs['OutputRate'] = 1.0
        fl_new['rel_vort'].attrs['ValidRange'] = [0.0, 1000.0]

    # add a special attribute if data is corrected explaining how this was done!
    if correct_data:
        disclaimer = 'These data have been limited to P-3 heights between 2.5 km and 3.5 km.\nThis allows for accurate comparisons.\nUnphysical outlier data have also been removed.'
        fl_new.attrs[ 'Note'] = disclaimer

    return fl_new


# just like the function above, but it trims down
def save_one_fl_processed( inputdata, limits, name):
    # add trimmed fields to this dataset!
    dvs = {}
    # make a field list for adding things in a separate loop
    field_list = []
    # fieldval is the NAME of the dataset, not the actual data...
    for fieldi, fieldval in enumerate( inputdata.keys()):
        # better case: correct for inconsistencies in the flight level data! will improve stats later
        field = inputdata[ fieldval] [ limits[0] : limits[1] ].values
        # field_list.append( field)
        # for fieldi, fieldval in enumerate( field_list):
        # add this field to the dataframe!
        dvs.update( { fieldval: ('time', field) })


    # next step: create the dataset and update units for existing variables!
    fl_new = xr.Dataset( data_vars=dvs, coords={'time': inputdata.time[ limits[0]:limits[1]].values })
    for fieldi, fieldval in enumerate( fl_new.variables):
        if fieldval in inputdata.variables:
            attribs = inputdata[ fieldval].attrs
            # add units, etc metadata present in the original flight level dataset!
            if 'units' in attribs:
                fl_new[fieldval].attrs['units'] = inputdata[ fieldval].units
            if 'Description' in attribs:
                fl_new[fieldval].attrs['Description'] = inputdata[ fieldval].Description
            if 'SampleRate' in attribs:
                fl_new[fieldval].attrs['SampleRate'] = inputdata[ fieldval].SampleRate
            if 'OutputRate' in attribs:
                fl_new[fieldval].attrs['OutputRate'] = inputdata[ fieldval].OutputRate
            if 'ValidRange' in attribs:
                fl_new[fieldval].attrs['ValidRange'] = inputdata[ fieldval].ValidRange
            if 'Min' in attribs:
                fl_new[fieldval].attrs['Min'] = inputdata[ fieldval].Min
            if 'Max' in attribs:
                fl_new[fieldval].attrs['Max'] = inputdata[ fieldval].Max


    # add a new attribute explaining how this dataset has been edited
    author = 'Created by NOAA HRD. Link to data: https://www.aoml.noaa.gov/ftp/hrd/data/flightlevel/'
    fl_new.attrs['Author'] = author
    disclaimer = 'Edited by: Ethan Murray (etmu9498@colorado.edu)'
    fl_new.attrs[ 'Editor'] = disclaimer
    disclaimer = 'This Dataset is based off the dataset ' + name + ". Please see that file for the original data."
    fl_new.attrs[ 'Attribution'] = disclaimer
    disclaimer = 'These data have been limited to P-3 heights between 2.5 km and 3.5 km.\nThis allows for accurate comparisons.\nUnphysical outlier data have also been removed.'
    fl_new.attrs[ 'Note'] = disclaimer

    return fl_new
