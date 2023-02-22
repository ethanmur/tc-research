# import packages and plotting scripts
# import packages and plotting scripts

import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
from datetime import datetime
import warnings
import pandas as pd
import cartopy.crs as ccrs
import shapely.geometry as sgeom
import matplotlib.patches as mpatches
from scipy.signal import find_peaks

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots


def load_in_situ( flight_path, flight_number, sample_step_size= 1):
    # make in situ data more manageable to work by using pandas, then convert to xarray
    # load in situ data

    os.chdir( flight_path)
    in_situ_data = pd.read_csv( flight_number, header=None, on_bad_lines='skip')

    # in_situ_data

    # make proper row the df header (some rows have text descriptors, for some reason this varies
    # between data sets)
    for col_ind in range( len( in_situ_data.columns)):
        if in_situ_data.iloc[ col_ind][1] == 'TIME':
            in_situ_data.columns = in_situ_data.iloc[ col_ind]
            break

    # removing 1203 rows that just repeat the header keys, and 6 rows with empty lat lon
    in_situ_data = in_situ_data[ in_situ_data['TIME'] != 'TIME']
    in_situ_data = in_situ_data[ ~ pd.isna( in_situ_data['LATref']) ]
    # removing 4 columns that are labeled as 'none'
    in_situ_data.drop( 'none', inplace=True, axis=1)

    # trim out every sample_step_size element for a new, smaller dataset. Easier / faster to work with
    # a step size of 1 won't change the data at all
    in_situ_data_trim = in_situ_data.iloc[ ::sample_step_size, :]
    # reset indices so they're nice and pretty!
    in_situ_data_trim = in_situ_data_trim.reset_index( drop=True)
    # in_situ_data_trim

    # adding datetime and just time formatted columns to pandas dataframe
    in_situ_data_trim['dt'] = pd.to_datetime( in_situ_data_trim['TIME'])
    in_situ_data_trim['time'] = [dt_object.time() for dt_object in in_situ_data_trim.dt ]

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

    # pull pitch and roll values out of the dataset, and save them as floats
    # I think they were originally strings
    pitch = xr_in_situ.PITCHref.values
    pitch = [ float( line) for line in pitch]

    roll = xr_in_situ.ROLLref.values
    roll = [ float( line) for line in roll]

    # newlon = [ float( value) for value in in_situ.LONref]
    # newlat = [ float( value) for value in in_situ.LATref]


    # add new time variables to dataset
    # older method: xr_in_situ.assign( str_time=str_time, float_time=float_time)
    xr_in_situ['float_time'] = ( ['index'], float_time)
    xr_in_situ['str_time'] = ( ['index'], str_time)
    xr_in_situ['rollval'] = ( ['index'], roll)
    xr_in_situ['pitchval'] = ( ['index'], pitch)

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

    # return the processed data for graphing use
    return xr_in_situ
