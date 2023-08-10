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
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023


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
    yearlist, filelist = helper_fns_winter2023.get_crl_datasets( tc, crl_data_root)

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
    print(crl_name)
    
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
        wv_2d = wv_2d.where( wv_2d.values < 30.)
        wv_2d = wv_2d.where( wv_2d.values >= 0.)
        lsr_2d = crl_data.LSR.where( crl_data.LSR.values < 10)
        lsr_2d = lsr_2d.where( lsr_2d.values > .1)

        # calculate a few more quantities that are a bit more involved!
        # pressure = find_crl_distance_rmws.find_p( crl_data, fl_data)
        # rh = find_crl_distance_rmws.find_rh( crl_data, fl_data)
        # temp_anoms = find_crl_distance_rmws.find_temp_anom()

        # use a helper fn to interpolate the data with the new height lims!
        vars = [ T_2d, power_2d, wv_2d, lsr_2d]
        var_names = ["T", "P_ch1", "WVMR", "LSR"]

        # special case for 9/8/22 earl data: fill in a roughly hour long time gap!
        if crl_name == "P3_20220908H1_095405-141758.cdf":
            print("Special Case for TC Earl: fill time gaps with nans")

            # find the time gap here
            starti = 0
            endi = 0
            time = crl_data.time.values
            for timei in range(len(time) - 1):
                if time[timei+1] - time[timei] > .1:
                    starti = timei
                    endi = timei+1
            # find an updated time matrix
            step = 2 / 3600 # 2 second timestep for crl data
            missingtimes = np.arange(time[starti] + step, time[endi], step) # increment the start by step to not duplicate that value
            crltimes = np.concatenate(( np.concatenate((time[0:starti], missingtimes), axis=0), time[endi:len(time)-1]), axis=0)

            '''
            print('start end times')
            print(starti)
            print(endi)
            print('crl times')
            print(len(crltimes))
            print("missing times info:")
            print(len(missingtimes))           
            print(missingtimes)
            print("Length of corrected time index")
            print(len(crltimes))
            print(crltimes)
            print(starti)
            print(endi)
            '''

            # correct the p3 heights on this step to match extended time axis
            p3_heights = find_crl_distance_rmws.find_interp( crl_data, fl_data, fl_data['HT.d'], correct_case=len(crltimes))

            # print('p3 heights')
            # print(len(p3_heights))

            # do this for each of the 4 matrices
            for fieldi, fieldval in enumerate( vars):
                nanarray = np.empty( (len(missingtimes), np.shape(fieldval)[1]))
                nanarray[:] = np.nan

                # now, slice the nan array into the current array
                array1, array2 = fieldval[0:starti, :], fieldval[endi:np.shape( fieldval)[0]-1, :]
                vars[fieldi] = np.concatenate(( np.concatenate((array1, nanarray), axis=0), array2), axis=0)
                totalarray = np.concatenate(( np.concatenate((array1, nanarray), axis=0), array2), axis=0)

                '''
                # print(len(crltimes))
                # print(np.shape(vars[fieldi]))
                print("Shape of 1st and second halves of temperature array:")
                print(np.shape(array1))
                print(np.shape(array2))
                print("shape of corrected array:")
                print(np.shape(vars[fieldi]))
                '''

        # otherwise, use default values
        else:
            crltimes = crl_data.time.values

        vars_interpolated = []
        for vari, varval in enumerate( vars):
            newh, var_2d = find_crl_distance_rmws.interp_data( varval, crl_data.H, p3_heights, crltimes,  year=yearval, crl_name=crl_name)
            # save
            if vari == 0:
                vars_interpolated.append( newh)
            vars_interpolated.append( var_2d)
            print( var_names[ vari] + " Interpolated")
            
    # add downscaled flight level data like wind speeds, etc!
    if add_dist_coords['fl_fields']:

        fl_field_list = [ 'WS.d', 'UWZ.d', 'MR.d', 'TA.d', 'PSURF.d']
        fl_field_interp = []
        for fieldi, fieldval in enumerate( fl_field_list):

            if crl_name == "P3_20220908H1_095405-141758.cdf":
                fl_field_interp.append( find_crl_distance_rmws.find_interp( crl_data, fl_data, fl_data[ fieldval], correct_case=len(crltimes)) )
            else:
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
    # special case for TC Julia: correct the time axis 
    # (remove the one problem data point identified in "tests\2023-05-18 tc julia axis correction.ipynb")
    time = crl_data.time.values
    if yearval == '2022' and crl_name[0:9] == "P3_202210":
        print("Julia Case!")
        # fill the outlier data point location with a nan
        time[np.where(time > 28.)[0]] = np.nan
        # interpolate the nan position
        nans, x = np.isnan(time), lambda z: z.nonzero()[0]
        time[nans] = np.interp(x(nans), x(~nans), time[~nans])

    # add height coords here
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

    # correct lat, lot, and radial distances here
    if crl_name == "P3_20220908H1_095405-141758.cdf":
        print("Special Case for TC Earl: correct other values")
        lonorig, latorig = crl_data.Lon.values, crl_data.Lat.values
        nanlist = np.empty( len(missingtimes))
        nanlist[:] = np.nan
        lons = np.concatenate(( np.concatenate((lonorig[0:starti], nanlist), axis=0), lonorig[endi:len(time)-1]), axis=0)
        lats = np.concatenate(( np.concatenate((latorig[0:starti], nanlist), axis=0), latorig[endi:len(time)-1]), axis=0)
        dvs={'ProductionDateTime':( crl_data.ProductionDateTime.values),'VersionID':( crl_data.VersionID.values),
                     'Lon':( 'time', lons), 'Lat':( 'time', lats) }

        # crldist = np.concatenate(( np.concatenate((crldist[0:starti], nanlist), axis=0), crldist[endi:len(time)-1]), axis=0)
        # crlrmw = np.concatenate(( np.concatenate((crlrmw[0:starti], nanlist), axis=0), crlrmw[endi:len(time)-1]), axis=0)
        crldist, crlrmw = find_crl_distance_rmws.find_rmws( crl_data, fl_data, correct_case=len(crltimes))

    else:
        dvs={'ProductionDateTime':( crl_data.ProductionDateTime.values),'VersionID':( crl_data.VersionID.values),
        	     'Lon':( 'time', crl_data.Lon.values), 'Lat':( 'time', crl_data.Lat.values) }

    # if new matrices were created, add them here!
    if add_dist_coords['new_heights']:
        temp = vars_interpolated[ 1]
        power = vars_interpolated[ 2]
        wv = vars_interpolated[ 3]
        lsr = vars_interpolated[ 4]
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
        units += ( 'm/s', 'm/s', "g/kg", "Degrees C", "hPa" )
        description1 += ( "Flight Level Total Wind Speed", "Flight Level Vertical Velocity", "Flight level water vapor mixing ratio",
                "Flight level temperature", "Flight level surface pressure")
        description2 += ( "", "", "", "", "")
        description3 += ( "", "", "", "", "")
        dvs.update( { 'wind_speed': ('time', fl_field_interp[0]), 'w': ('time', fl_field_interp[1]), 'fl_wv': ( 'time', fl_field_interp[2]),
        'fl_T': ( 'time', fl_field_interp[3]), 'fl_psurf': ('time', fl_field_interp[4])})


    # next step: create the dataset and update units!
    crl_new = xr.Dataset( data_vars=dvs, coords={'time': crltimes, 'height': height})

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
