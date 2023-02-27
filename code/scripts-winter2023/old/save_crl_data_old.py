# import...
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings
import sys
os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/crl-data-processing")
import find_crl_distance_rmws



# this function is the main loop to save multiple new CRL datasets.
# inputs:
############
# tc: can either be 'all' to do this for all 2021 and 2022 data, a year string like
#     '2021' or '2022', or a dictionary with years as the key and files to look at in
#     a list as the values. The last option is the most useful for smaller tests!
# add_dist_coords: a dictionary of possible distance coordinates to add / other optional
#     inputs. All the inputs default to False (kinda useless lol but prevents errors!)
def save_tcs( tc='all', add_dist_coords={'new_heights': False, 'fl_fields': False, 'rmw': False}):

    crl_data_root = "/Users/etmu9498/research/data/CRL_data/"

    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']

        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( make_plots.load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( make_plots.load_flight_level( crl_data_root + '2022', print_files=False) )

    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2022', print_files=False)]

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
            crl_data = save_one_crl( yearval, fileval, add_dist_coords=add_dist_coords)

            # make a name for the dataset
            filename = filelist[ yeari][ filei] [ : -18] + "_processed.nc"
            # save the newly created crl dataset
            crl_data.to_netcdf('/Users/etmu9498/research/data/crl-all-data-processed/' + yearval + "/" + filename)

            print( "New CRL File Created and Saved: " +  yearval + "/" + filename)
    return



def save_one_crl( yearval, crl_name, add_dist_coords={'new_heights': False, 'fl_fields': False, 'rmw': False}):

    # open up the crl file
    crl_data_root = "/Users/etmu9498/research/data/CRL_data/"
    crl_path = crl_data_root + yearval
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # make a copy of the crl data: this will eventually become the new, full dataset!
    crl_new = crl_data.copy()
    # add a new attribute explaining how this dataset has been edited
    disclaimer = 'This Dataset is based off the dataset ' + crl_name + ". Please see that file for the original data."
    crl_new.attrs[ 'global_att4'] = disclaimer


    # go through every key for add_dist_coords. If any of them are true, we'll need to
    # load the crl's corresponding flight level dataset!
    load_data = False
    for keyval in add_dist_coords.values():
        if keyval:
            load_data = True
    if load_data:
        fl_data = find_crl_distance_rmws.find_matching_fl( yearval, crl_name)
        print( "Corresponding flight level dataset loaded")

    # set time as the only true coordinate!
    crl_new = crl_new.set_coords( 'time')

    if add_dist_coords['rmw']:

        # downsample the radial distance and rmw coordinates found in the flight level dataset
        # to this new dataset!
        crldist, crlrmw = find_crl_distance_rmws.find_rmws( crl_data, fl_data)

        print( "In Situ radial distance and rmw axes created")

        # save the in situ distance array in a separate variable
        crl_new = crl_new.assign( {'center_dist': crldist }, dims=['time'])
        # crl_new[ 'center_dist'] = crldist
        crl_new.center_dist.attrs['long_name'] = 'center_dist'
        crl_new.center_dist.attrs['units'] = 'km'
        crl_new.center_dist.attrs['description'] = 'Distance from the center of the TC.'
        crl_new.center_dist.attrs['description2'] = 'Created using P-3 locations and NOAA official TC tracks.'
        crl_new.center_dist.attrs['description3'] = 'nan values represent missing or non-overlapping P-3 and / or TC track locations.'

        crl_new = crl_new.assign( {'rmw': crlrmw })
        # crl_new[ 'rmw'] = crlrmw
        crl_new.rmw.attrs['long_name'] = 'rmw'
        crl_new.rmw.attrs['units'] = 'unitless'
        crl_new.rmw.attrs['description'] = 'The Radius of Maximum Winds for every eye pass.'
        crl_new.rmw.attrs['description2'] = 'Created using maximum eyewall wind speeds and the radial distances described above.'
        crl_new.rmw.attrs['description3'] = 'nan values represent missing or non-overlapping P-3 and / or TC track locations, or TCs without peak wind speeds.'

        # reset radial distance and rmws as variables, not coordinates!
        # crl_new = crl_new.reset_coords( [ 'rmw', 'center_dist'])


    # find the p-3 heights using a helper function, then use another helper function to
    # interpolate the data!
    if add_dist_coords['new_heights']:

        # make a copy of the crl dataset: variables below will be dropped to make room for updated data!
        # important: need to drop "H" before adding a new one! or else xarray won't be happy
        crl_new = crl_new.drop_vars( ["H", "T", "P_ch1", "WVMR", 'LSR'])

        # first, find the p-3 height values using interpolation in a helper function
        p3_heights = find_crl_distance_rmws.find_p3height( crl_data, fl_data)
        maxh = np.nanmax( p3_heights)
        # crl_new = crl_new.assign_coords( {'H_max': maxh })
        crl_new["P3_height"] = p3_heights
        crl_new.rmw.attrs['units'] = 'm from surface'
        crl_new.rmw.attrs['description'] = 'The P-3 height at any given time taken from flight level data.'

        # fix matrices before running the helper function to account for height
        # based on code in the make_plots function
        ## important!!!
        # these are old settings: feel free to change!
        T_2d = crl_data.T.where( crl_data.T.values < 50)

        power_2d = 10 * np.log10( crl_data.P_ch1 )
        power_2d = power_2d.where( power_2d.values > -40) # -30

        wv_2d = crl_data.WVMR.where( crl_data.WVMR.values != 0)
        wv_2d = wv_2d.where( wv_2d.values < 30)

        lsr_2d = crl_data.LSR.where( crl_data.LSR.values < 10)
        lsr_2d = lsr_2d.where( lsr_2d.values > .1)


        # calculate a few more quantities that are a bit more involved!
        # pressure = find_crl_distance_rmws.find_p( crl_data, fl_data)
        # rh = find_crl_distance_rmws.find_rh( crl_data, fl_data)
        # temp_anoms = find_crl_distance_rmws.find_temp_anom()

        # use a helper fn to interpolate the data with the new height lims!
        vars = [ T_2d, power_2d, wv_2d, lsr_2d]
        var_names = ["T", "P_ch1", "WVMR", "LSR"]

        for vari, varval in enumerate( vars):
            newh, var_2d = find_crl_distance_rmws.interp_data( varval, crl_data.H, p3_heights)
            crl_new = crl_new.assign( { var_names[vari] : xr.DataArray( var_2d)})

        crl_new["H"] = - newh * 1000
        crl_new.rmw.attrs['units'] = 'm from surface'
        crl_new.rmw.attrs['description'] = 'The interpolated height array for all CRL measurement matrices.'

    return crl_new
