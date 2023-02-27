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


    # go through every key for add_dist_coords. If any of them are true, we'll need to
    # load the crl's corresponding flight level dataset!
    load_data = False
    for keyval in add_dist_coords.values():
        if keyval:
            load_data = True
    if load_data:
        fl_data = find_crl_distance_rmws.find_matching_fl( yearval, crl_name)
        print( "Corresponding flight level dataset loaded")


    #############
    ## step 1: use flags to find updated / new arrays and matrices!
    ##         don't add them yet tho... that's the next step!
    #############

    if add_dist_coords['rmw']:
        # downsample the radial distance and rmw coordinates found in the flight level dataset
        # to this new dataset!
        crldist, crlrmw = find_crl_distance_rmws.find_rmws( crl_data, fl_data)
        print( "In Situ radial distance and rmw axes found")

    # find the p-3 heights using a helper function, then use another helper function to
    # interpolate the data!
    if add_dist_coords['new_heights']:
        # first, find the p-3 height values using interpolation in a helper function
        p3_heights = find_crl_distance_rmws.find_interp( crl_data, fl_data, fl_data['HT.d'])

        # fix matrices before running the helper function to account for height
        # based on code in the make_plots function
        ## important!! these are old settings: feel free to change!
        T_2d = crl_data.T.where( crl_data.T.values < 50)
        warnings.filterwarnings("ignore")
        power_2d = 10 * np.log10( crl_data.P_ch1 )
        power_2d = power_2d.where( power_2d.values > -40)
        power_2d = power_2d.where( power_2d.values < 0)
        warnings.filterwarnings("default")
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

        vars_interpolated = []
        for vari, varval in enumerate( vars):
            newh, var_2d = find_crl_distance_rmws.interp_data( varval, crl_data.H, p3_heights)
            # save
            if vari == 0:
                vars_interpolated.append( newh)
            vars_interpolated.append( var_2d)

    # add downscaled flight level data like wind speeds, etc!
    if add_dist_coords['fl_fields']:

        fl_field_list = [ 'WS.d', 'UWZ.d', 'MR.d', 'TA.d' ]
        fl_field_interp = []
        for fieldi, fieldval in enumerate( fl_field_list):
            fl_field_interp.append( find_crl_distance_rmws.find_interp( crl_data, fl_data, fl_data[ fieldval]) )




    ###########
    ## step 2: create a new xarray dataset and save the calculated values!
    ##         this is done after caclulateing all the variables because one of the
    ##         coords (height) changes if add_dist_coords is true or false
    ###########

    # use these lists to recursively go through each case
    units = []
    description1 = []
    description2 = [ ]
    description3 = [ ]
    #making time and height coordinates
    time = crl_data.time.values
    if add_dist_coords['new_heights']:
        height = - 1000 * vars_interpolated[ 0]
    else:
        height = - 1000 * crl_data.H
    sv_dims=['time', 'height']


    # add the variables below no matter what!
    units += ( "", "", "Degrees West", "Degrees North")
    description1 += ( "", "", "", "")
    description2 += ( "", "", "", "")
    description3 += ( "", "", "", "")
    dvs={'ProductionDateTime':( crl_data.ProductionDateTime.values),'VersionID':( crl_data.VersionID.values),
    	     'Lon':( 'time', crl_data.Lon.values), 'Lat':( 'time', crl_data.Lat.values) }

    # if new matrices were created, add them here!
    if add_dist_coords['new_heights']:
        temp = vars_interpolated[ 1]
        power = vars_interpolated[ 2]
        wv = vars_interpolated[ 3]
        lsr = vars_interpolated[ 4]

        #print( len( time))
        #print( len( height))
        #print( np.shape( temp))
        #print( type( temp))
        #print( temp)

        # add metadata
        units += ( "Degrees C", "dBz", "g/kg", "unitless", 'm')
        description1 += ( "Temperature", "Returned channel 1 power", "Water vapor mixing ratio", "Light scattering ratio", "P-3 height above surface")
        description2 += ( "", "", "", "", "")
        description3 += ( "", "", "", "", "")
        dvs.update( { 'T': (sv_dims, temp), 'P_ch1': (sv_dims, power),
                'WVMR': (sv_dims, wv), 'LSR': (sv_dims, lsr) ,
                'p3_height': ( 'time', p3_heights) })

    # otherwise, add the old, original matrices!
    else:
        units += ( "Degrees C", "dBz", "g/kg", "unitless")
        description1 += ( "Temperature", "Returned channel 1 power", "Water vapor mixing ratio", "Light scattering ratio")
        description2 += ( "", "", "", "")
        description3 += ( "", "", "", "")
        dvs.update( { 'T': (sv_dims, crl_data.T), 'P_ch1': (sv_dims, crl_data.P_ch1),
        'WVMR': (sv_dims, crl_data.WVMR), 'LSR': (sv_dims, crl_data.LSR) })

    # if center distances and rmws were created, add them in this step!
    if add_dist_coords['rmw']:
        units += ( "Km", "unitless")
        description1 += ( 'Distance from the center of the TC.', 'The Radius of Maximum Winds for every eye pass.')
        description2 += ( 'Created using P-3 locations and NOAA official TC tracks.',
                "Created using maximum eyewall wind speeds and the radial distances described above.")
        description3 += ( 'nan values represent missing or non-overlapping P-3 and / or TC track locations.',
                'nan values represent missing or non-overlapping P-3 and / or TC track locations, or TCs without peak wind speeds.')
        dvs.update( { 'center_dist': ( 'time', crldist), 'rmw': ( 'time', crlrmw) })


    # add the new downscaled flight level data!
    if add_dist_coords['fl_fields']:
        units += ( 'm/s', 'm/s', "g/kg", "Degrees C" )
        description1 += ( "Flight Level Total Wind Speed", "Flight Level Vertical Velocity", "Flight level water vapor mixing ratio",
                "Flight level temperature")
        description2 += ( "", "", "", "")
        description3 += ( "", "", "", "")
        dvs.update( { 'wind_speed': ('time', fl_field_interp[0]), 'w': ('time', fl_field_interp[1]), 'fl_wv': ( 'time', fl_field_interp[2]),
        'fl_T': ( 'time', fl_field_interp[3]) })


    # next step: create the dataset and update units!
    crl_new = xr.Dataset( data_vars=dvs, coords={'time': time, 'height': height})

    # add explanations to global attributes
    # add a new attribute explaining how this dataset has been edited
    disclaimer = 'This Dataset is based off the dataset ' + crl_name + ". Please see that file for the original data."
    crl_new.attrs[ 'global_att1'] = crl_data.attrs[ 'global_att']
    crl_new.attrs[ 'global_att2'] = crl_data.attrs[ 'global_att2']
    crl_new.attrs[ 'global_att3'] = crl_data.attrs[ 'global_att3']
    crl_new.attrs[ 'global_att4'] = disclaimer


    # add units and explanations for each variable
    for keyi, keyval in enumerate( dvs.keys()):
        crl_new[ keyval].attrs['units'] = units[ keyi]
        if description1[ keyi] != '':
            crl_new[ keyval].attrs['description'] = description1[ keyi]
        if description2[ keyi] != '':
            crl_new[ keyval].attrs['description2'] = description2[ keyi]
        if description3[ keyi] != '':
            crl_new[ keyval].attrs['description3'] = description3[ keyi]

    # add coordinate metadata and units!
    crl_new[ 'time'].attrs['units'] = 'hours (UTC, decimal)'
    crl_new[ 'height'].attrs['units'] = 'm from surface'
    crl_new[ 'height'].attrs['description'] = 'linear height coordinate for 2D data'

    return crl_new
