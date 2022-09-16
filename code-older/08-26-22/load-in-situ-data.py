# Author: Ethan Murray
# Date: 9/8/22

# methods used to load and process in situ data

import pandas as pd
import os
import xarray as xr

flight_path = "/Path/to/data"
flight_name = '20210816H1_iwg1.txt'

# load data-
os.chdir( flight_path)
in_situ_data = pd.read_csv( flight_name, header=None)

# some rows have text descriptors instead of data, for some reason their location varies between data sets)
# make one of these rows the header, and get rid of the rest
for row_ind in range( len( in_situ_data.columns)):

    # this finds the text descriptor rows, and makes one of them the header
    if in_situ_data.iloc[ row_ind][1] == 'TIME':
        in_situ_data.columns = in_situ_data.iloc[ row_ind]
        break
# removing 1203 rows that just repeat the header keys, and 6 rows with empty lat lon
in_situ_data = in_situ_data[ in_situ_data['TIME'] != 'TIME']
in_situ_data = in_situ_data[ ~ pd.isna( in_situ_data['LATref']) ]
# removing 4 columns that are labeled as 'none'
in_situ_data.drop( 'none', inplace=True, axis=1)


# trim out every sample_step_size element for a new, smaller dataset.
# It's sometimes easier and faster to work with a smaller dataset, but a step size of 1 won't change the data at all
sample_step_size = 1
in_situ_data_trim = in_situ_data.iloc[ ::sample_step_size, :]
# reset indices
in_situ_data_trim = in_situ_data_trim.reset_index( drop=True)

# adding datetime and just time formatted columns to pandas dataframe
# these different time formats are helpful for plotting and calculating timesteps
in_situ_data_trim['dt'] = pd.to_datetime( in_situ_data_trim['TIME'])
in_situ_data_trim['time'] = [dt_object.time() for dt_object in in_situ_data_trim.dt ]

# convert from pandas to xarray
xr_in_situ =  pd.DataFrame.to_xarray( in_situ_data_trim)

# get data out of xarray and put it in a useable format
time = xr_in_situ.time.values
# create string and decimal time variables for easier plotting
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
# you can't manipulate the data at all if they're saved as strings
pitch = xr_in_situ.PITCHref.values
pitch = [ float( line) for line in pitch]

roll = xr_in_situ.ROLLref.values
roll = [ float( line) for line in roll]

# add new time variables to dataset
# older method: xr_in_situ.assign( str_time=str_time, float_time=float_time)
xr_in_situ['float_time'] = ( ['index'], float_time)
xr_in_situ['str_time'] = ( ['index'], str_time)
xr_in_situ['rollval'] = ( ['index'], roll)
xr_in_situ['pitchval'] = ( ['index'], pitch)

# this line turns all coordinates into variables
# not really necessary, but I like having everything as a variable
xr_in_situ.reset_coords()
