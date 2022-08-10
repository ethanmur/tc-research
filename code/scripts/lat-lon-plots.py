import os
import xarray as xr
import matplotlib.pyplot as plt

os.chdir( "/Users/etmu9498/research/code/scripts/in-situ-plots")
import load_in_situ_data

def plot_lat_lon( crl_path, crl_name, flight_data_path, flight_name, i1, i2, xaxis, variable_list, xlim='none', title=False):

    # load and process the data
    xr_in_situ = load_in_situ_data.load_in_situ( flight_data_path, flight_name, sample_step_size=10)

    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]

    plt.figure( figsize=(8, 8))

    start_ind = crl_data.time[ i1]
    end_ind = crl_data.time[ i2]
    plt.axvline( x= start_ind, c='g')
    plt.axvline( x= end_ind, c='g')
    if xaxis == "lon":
        xaxis_data = lon
        plt.plot( float_time, lon )
        xlabel = 'Longitude (Degrees)'
    elif xaxis == "lat":
        xaxis_data = lat
        plt.plot( float_time, lat)
        xlabel = 'Latitude (Degrees)'

    plt.xlabel( "Time (UTC, Hours)")
    plt.ylabel( xlabel)
    plt.title( title)
    plt.grid('on')
