import os
import pandas as pd
import xarray as xr
import warnings

os.chdir( "/Users/etmu9498/research/code/scripts")
import make_plots

def save_as_netcdf():

    # setup
    txt_path = "/Users/etmu9498/research/data/in-situ"
    os.chdir( txt_path)
    txt_list = make_plots.load_flight_level( txt_path, print_files=False)

    # do this for every in situ data file
    for dataset in range( len( txt_list)):
        warnings.filterwarnings("ignore")

        in_situ_data = pd.read_csv( txt_list[dataset], header=None, on_bad_lines='skip')

        print( 'Pandas dataset loaded for file ' + str( dataset))

        # make proper row the df header (some rows have text descriptors, for some reason this varies
        # between data sets)
        for col_ind in range( len( in_situ_data.columns)):
            # print( col_ind)
            if in_situ_data.iloc[ col_ind][1] == 'TIME':
                in_situ_data.columns = in_situ_data.iloc[ col_ind]
                break

        # removing 1203 rows that just repeat the header keys, and 6 rows with empty lat lon
        in_situ_data = in_situ_data[ in_situ_data['TIME'] != 'TIME']
        in_situ_data = in_situ_data[ ~ pd.isna( in_situ_data['LATref']) ]
        # removing 4 columns that are labeled as 'none'
        in_situ_data.drop( 'none', inplace=True, axis=1)

        # new step!
        # cut pandas dataframe down to only the values that are of interest
        # this will hopefully save some time when saving and transitioning to xr?
        # keylist was taken from simple_flight_level_plot.plot()... shoud be all I need!
        keylist = [ 'TIME', 'WS.d', 'WD.d', 'UWZ.d', 'ASfmrRainRate.1', 'LATref', 'LONref', 'TAS.d', 'PSURF.d', 'HT.d', 'PITCHref', 'ROLLref']
        in_situ_data_trim = in_situ_data[ keylist]

        print( "Pandas dataset trimmed for file " + str( dataset))

        # adding datetime and just time formatted columns to pandas dataframe
        in_situ_data_trim['dt'] = pd.to_datetime( in_situ_data_trim['TIME'])
        in_situ_data_trim['time'] = [dt_object.time() for dt_object in in_situ_data_trim.dt ]

        # this snippet turns all relevant quantities from strings into floats for easier use later
        keylist2 = [ 'WS.d', 'WD.d', 'UWZ.d', 'ASfmrRainRate.1', 'LATref', 'LONref', 'TAS.d', 'PSURF.d', 'HT.d', 'PITCHref', 'ROLLref']
        for key in keylist2:

            in_situ_data_trim[ key] = in_situ_data_trim[ key].astype( float)
            # print( 'Key ' + key + ' saved')



        # convert from pandas to xarray
        xr_in_situ =  pd.DataFrame.to_xarray( in_situ_data_trim)

        # get data out of xarray and put it in a useable format
        time = xr_in_situ.time.values
        # also store the data in string and decimal formats for easier plotting
        str_time = [ ti.strftime('%H:%M:%S') for ti in time ]
        float_time = []
        for val in time:
            # account for wraparound times into next day
            if val.hour + val.minute / 60 + val.second / 3600 < 6.0:
                float_time.append( 24.0 + val.hour + val.minute / 60 + val.second / 3600)
            else:
                float_time.append( val.hour + val.minute / 60 + val.second / 3600)
        # xr_in_situ['float_time'] = ( ['index'], float_time)
        xr_in_situ['str_time'] = ( ['index'], str_time)

        # this variable has an ambiguous 'object' data type, which causes problems
        # when trying to save the dataset later... so it needs to be dropped here!
        xr_in_situ = xr_in_situ.drop( 'time')


        # change time from a variable to a coordinate
        # use float_time's values to avoid ambiguous data types
        xr_in_situ= xr_in_situ.assign_coords( {'time': float_time })
        xr_in_situ.time.attrs['long_name'] = 'time'
        xr_in_situ.time.attrs['units'] = 'Hours (UTC)'
        xr_in_situ.time.attrs['description'] = "Time of measurement, from the P-3's internal clock"

        # this line ...
        xr_in_situ.reset_coords()

        '''
        # this snippet turns all relevant quantities from strings into floats for easier use later
        keylist2 = [ 'WS.d', 'WD.d', 'UWZ.d', 'ASfmrRainRate.1', 'LATref', 'LONref', 'TAS.d', 'HT.d', 'PITCHref', 'ROLLref']
        for key in keylist2:

            floatkey = xr_in_situ[ key].astype( float)

            # floatkey = [ float( line) for line in xr_in_situ[key]]
            xr_in_situ[ key] = ('index', floatkey)
            print( 'Key ' + key + ' saved')




            key_temp = np.zeros( len( xr_in_situ[ key]))
            count = 0
            for line_ind in range( len( return_var)):
                if return_var[ line_ind] == '':
                    return_var_temp[line_ind] = np.nan
                    count += 1
                else:
                    return_var_temp[ line_ind] = float( return_var[ line_ind])
            return_var_temp.tolist()

        '''

        warnings.filterwarnings("default")

        # save the new datasets!
        filename = txt_list[dataset][0:-3] + 'nc'
        xr_in_situ.to_netcdf('/Users/etmu9498/research/data/in-situ-nc/' + filename)
        print( "In Situ File Created and Saved: " + filename)
