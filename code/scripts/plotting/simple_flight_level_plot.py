import os
import warnings
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
os.chdir("/Users/etmu9498/research/code/scripts")
import tc_metadata

def plot(crl_path, crl_name, tcdata, counter, xaxis='dist', plot_number=313):
    warnings.filterwarnings("ignore")

    in_situ_path = tcdata[ 'new_flight_data_path']
    in_situ_name = tc_metadata.choose_new_in_situ_name( tcdata['tc_name'], counter)

    # load in situ data
    os.chdir( in_situ_path)
    xr_in_situ = xr.open_dataset( in_situ_name)
    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # crop in situ data
    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time
    time1 = crl_data.time[ 0].values
    time2 = crl_data.time[ -1].values
    # find the in situ times nearest the crl times
    idx1 = (np.abs(float_time - time1)).argmin()
    idx2 = (np.abs(float_time - time2)).argmin()


    keyList = [ 'WS.d', 'WD.d', 'UWZ.d', 'ASfmrRainRate.1', 'LATref', 'LONref', 'TAS.d', 'HT.d', 'distance']
    # make an empty dict that will be filled soon!
    datatrim = {key: None for key in keyList}
    # fill dict with keys provided after editing data!
    for key in keyList:
        temp_var = empty_str_helper( xr_in_situ[key].values)
        temp_var = temp_var[ idx1.values : idx2.values]
        datatrim[ key] = temp_var



    # get the proper axis
    if xaxis == "time":
        xaxis_data = float_time
        xlabel = 'Time (UTC, Hours)'
    elif xaxis == "lon":
        xaxis_data = datatrim['LONref']
        xlabel = 'Longitude (Degrees)'
    elif xaxis == "lat":
        xaxis_data = datatrim['LATref']
        xlabel = 'Latitude (Degrees)'
    elif xaxis == 'dist':
        xaxis_data = datatrim['distance']
    else:
        print( "Please choose 'lon', 'lat', 'time', or 'dist' as a valid xaxis!")

    # make sure the x axis data is a numpy array for proper plotting
    xaxis_data = np.array(xaxis_data)

    fig = plt.gcf()
    ax1 = fig.add_subplot( plot_number)

    ax1.plot( xaxis_data, datatrim['WS.d'], c='c', label='Tangential Wind Speed (m/s)')
    ax1.set_ylabel('Total Wind Speed (m/s)', c='c')
    ax1.xaxis.grid( )
    ax1.yaxis.grid( )

    ax2 = ax1.twinx()
    # ax2.plot( xaxis_data, datatrim['ASfmrRainRate.1'], c='b', label=' SFMR Rain Rate (mm/hr)')
    ax2.plot( xaxis_data, datatrim['UWZ.d'], c='y', label='Vertical Wind Speed ( m/s)')
    ax2.set_ylabel( 'Vertical Wind Speed (m/s', c='y')

    ax2.set_xlabel( 'Distance from TC Center (km)')

    # plt.legend(loc='upper left')
    # plt.grid('on')

    warnings.filterwarnings("ignore")

# turn all string values into floats
# this is a longer function to account for empty strings: turn those values into nans
def empty_str_helper( return_var):
    return_var_temp = np.zeros( len( return_var))
    count = 0
    for line_ind in range( len( return_var)):
        if return_var[ line_ind] == '':
            return_var_temp[line_ind] = np.nan
            count += 1
        else:
            return_var_temp[ line_ind] = float( return_var[ line_ind])
    return return_var_temp.tolist()

'''
    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time
    ws = [ float( line) for line in xr_in_situ["WS.d"].values] # convert strings to floats
    wd = [ float( line) for line in xr_in_situ["WD.d"].values]
    uwz = empty_str_helper( xr_in_situ["UWZ.d"].values)
    rr = empty_str_helper( xr_in_situ["ASfmrRainRate.1"].values)
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]
    tas = [float( line) for line in xr_in_situ["TAS.d"].values ]


    # use the indices for nearest time values to trim down lat and lon values
    # this prevents data overlap and makes plotting faster!
    lon = lon[ idx1.values : idx2.values]
    lat = lat[ idx1.values : idx2.values]
    ws = ws[ idx1.values : idx2.values]
    uwz = uwz[ idx1.values : idx2.values]
    rr = rr[ idx1.values : idx2.values]
    tas = tas[ idx1.values : idx2.values]
    dist = xr_in_situ.distance[ idx1.values : idx2.values]

'''
