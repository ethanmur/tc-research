## Last edited 5/27/22


import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import xarray as xr
import warnings
import datetime
from matplotlib.ticker import FormatStrFormatter


os.chdir(  "/Users/etmu9498/research/code/scripts/")
import helper_fns


# a helper script used in most plotting functions to determine the x axis scale variable

def x_axis_helper( data_path, data_file, index1, index2, xaxis):

    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    # choose x axis type
    if xaxis == 'lon':
        xaxis = crl_data.Lon[index1:index2]
        x_label = 'longitude (degrees)'
        xlims = [ crl_data.Lon[index1], crl_data.Lon[index2] ]
        return xaxis, x_label, xlims

    elif xaxis == 'lat':
        xaxis = crl_data.Lat[index1:index2]
        x_label = 'latitude (degrees)'
        xlims = [ crl_data.Lat[index1], crl_data.Lat[index2] ]
        return xaxis, x_label, xlims

    elif xaxis == 'time':
        xaxis = crl_data.time[index1:index2]
        x_label = 'Time (UTC)'
        xlims = [ crl_data.time[index1], crl_data.time[index2] ]
        return xaxis, x_label, xlims

    elif xaxis == 'dist':
        max_lat = np.max( crl_data.Lat[index1:index2] )
        min_lat = np.min( crl_data.Lat[index1:index2] )
        # using the delta lat formula found on the wikipedia page, assuming lat = 15 degrees
        scale = 110.649

        # make the 0 km mark the center of the tc!
        xaxis = np.linspace( 0, scale * ( max_lat - min_lat), index2 - index1)
        xaxis = xaxis - np.max( xaxis) / 2
        x_label = 'Distance (km)'
        xlims = [ np.min( xaxis), np.max( xaxis) ]
        return xaxis, x_label, xlims

    else:
        print("Error: Please Choose 'lat', 'lon', 'time', or 'distance' for the x axis")
        return


def plot_T(data_path, data_file, index1, index2, xaxis_name):

    warnings.filterwarnings("ignore")

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    # pick plot color
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()

    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate temperature
    temp = crl_data.T[index1:index2, :].where( crl_data.T[index1:index2, :].values < 50).transpose()

    # plot and make things pretty
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.pcolormesh( xaxis, - crl_data.H, temp, cmap = color_map )
    ax.set_ylabel( 'Height (km)')
    ax.set_xlabel( x_label)
    ax.set_xlim( xlims )
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.grid( )
    ax.yaxis.grid( )

    map = mpl.cm.ScalarMappable(cmap=color_map, norm=mpl.colors.Normalize(vmin=0, vmax=40))
    cbar = fig.colorbar(mappable = map, label="Temperature ( C)")

    return fig, ax
    warnings.filterwarnings("default")




def plot_T_anomaly(data_path, data_file, index1, index2, xaxis_name):

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate temperature and temperature anomaly
    temp = crl_data.T[index1:index2, :].where( crl_data.T[index1:index2, :].values < 50).transpose()

    # temperature anomaly is found using a local average from our data! An array
    layer_avg_temp = np.nanmean( crl_data.T.where( crl_data.T.values < 50).transpose(), axis = 1)
    # make an array of average temperatures within the correct bounds. A matrix
    avg_temp_obs = np.ones_like( temp)
    for i in range( index2 - index1):
        avg_temp_obs[:, i] = avg_temp_obs[:, i] * layer_avg_temp
    temp_anomaly_obs = temp - avg_temp_obs

    # set 0 as the central point of the pcolormesh!
    divnorm = colors.TwoSlopeNorm(vmin=-5, vcenter=0, vmax=15)
    # plot things and make pretty
    plt.pcolormesh( xaxis, - crl_data.H, temp_anomaly_obs, cmap = "bwr", norm=divnorm )
    plt.colorbar(label="Temperature Anomaly ( C)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')



def plot_wvmr(data_path, data_file, index1, index2, xaxis_name):

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate wvmr
    step1 = crl_data.WVMR.where( crl_data.WVMR.values != 0)
    step2 = step1.where( step1.values < 20)
    crl_wvmr = step2[index1:index2, :].transpose()

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H, crl_wvmr, vmin = 0, vmax =20)
    plt.colorbar(label="WVMR ( g/kg)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')



def plot_lsr(data_path, data_file, index1, index2, xaxis_name):

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate lsr
    step1 = crl_data.LSR[index1:index2, :].where( crl_data.LSR[index1:index2].values < 10).transpose()
    crl_lsr = step1.where( step1.values > .1)

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H, crl_lsr)
    plt.colorbar(label="LSR ( unitless)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')



def plot_power_ch1(data_path, data_file, index1, index2, xaxis_name, cutoff=-30):
    warnings.filterwarnings("ignore")

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate power backscattered to channel 1
    step1 = 10 * np.log10( crl_data.P_ch1 )

    step2 = step1.where( step1.values > cutoff)
    crl_pch1 = step2[index1:index2, :].transpose()


    # plot and make things pretty
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.pcolormesh(  xaxis, - crl_data.H, crl_pch1, vmin = cutoff, vmax =-10)
    ax.set_ylabel( 'Height (km)')
    ax.set_xlabel( x_label)
    ax.set_xlim( xlims )
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.grid( )
    ax.yaxis.grid( )

    color_map = 'viridis'
    map = mpl.cm.ScalarMappable(cmap=color_map, norm=mpl.colors.Normalize(vmin=cutoff, vmax=-10))
    cbar = fig.colorbar(mappable = map, label="Backscattered Ch 1 power ( dBz)")

    warnings.filterwarnings("default")
    return fig, ax, xlims



def plot_rh(data_path, data_file, index1, index2, xaxis_name):

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate rh
    # defining constants
    epsilon = .622
    e_0 = 6.112 # hPa
    b = 17.67
    T_1 = 273.15 # K
    T_2 = 29.65 # K
    # calculate wvmr
    step1 = crl_data.WVMR.where( crl_data.WVMR.values != 0)
    step2 = step1.where( step1.values < 20)
    wvmr = step2[index1:index2, :].transpose()
    # calculate temperature (K)
    temp = crl_data.T[index1:index2, :].where( crl_data.T[index1:index2, :].values < 50).transpose() + 273
    # define pressure field using scale height
    # no negative sign in exp() because heights are already negative
    scale_ht = 7.5 # km, just an estimate
    pressure = 1013.3 * np.exp( crl_data.H / scale_ht) # hPa
    # find saturation vapor pressure
    e_s = e_0 * np.exp(  ( b * ( temp - T_1) ) / (temp - T_2) )
    # find relative humidity!
    saturation_wvmr = 1000 * (epsilon * e_s ) / ( pressure - e_s) # multiply by 1000 to get to g/kg like wvmr (above)
    crl_rh = 100 * wvmr / saturation_wvmr
    crl_rh = crl_rh.where( crl_rh.values <= 100)
    crl_rh = crl_rh.where( crl_rh.values >= 0)

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H, crl_rh)
    plt.colorbar(label="RH ( %)")
    plt.ylabel( 'Height (km)')

    # unique fix for rh code axes, they need to be flipped for some reason to get plots to match :/
    if xaxis_name == 'lon' or xaxis_name == 'lat':
        plt.gca().invert_xaxis()

    # plt.xlabel( x_label)
    # plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')



def plot_sondes(crl_path, crl_file_name, sonde_path, variable_to_plot):
    # make a new height array for the rh measurements
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_file_name)

    os.chdir( sonde_path)
    # figure out how many files are in the selected sonde folder
    file_names = []
    for (dirpath, dirnames, file) in os.walk( sonde_path):
        file_names.extend(file)
        break

    # plot each dropsonde profile
    for ind in range( len( file_names)):
        sonde = xr.open_dataset( file_names[ ind])
        # look at the initial time for each sonde
        # do some horrible conversions to get result into a float, there's probably a better solution
        start_string_time = sonde.time[0].values.astype(str)[11:22]
        start_time = float( start_string_time[0:2]) + float( start_string_time[3:5]) / 60 + float( start_string_time[6:8]) / 3600
        end_string_time = sonde.time[-1].values.astype(str)[11:22]
        end_time = float( end_string_time[0:2]) + float( end_string_time[3:5]) / 60 + float( end_string_time[6:8]) / 3600
        mid_time = ( start_time + end_time ) / 2

        if mid_time < 3.0:
            mid_time = mid_time + 24

        # get temperatures
        Z = sonde[variable_to_plot].values[~np.isnan( sonde[variable_to_plot].values)]
        heights = np.linspace( np.nanmin( -crl_data.H.values), np.nanmax( -crl_data.H.values), len( Z))
        plt.pcolormesh( [mid_time - .00125, mid_time + .00125], heights, np.matrix( [Z, Z] ).transpose(),
                       vmin = 0, vmax =20 )




# the function I'm using to plot tdr reflectivity
def plot_tdr( tdr_path, inbound_name, outbound_name, xaxis):

    warnings.filterwarnings("ignore")

    # get data
    os.chdir( tdr_path)
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # choose x axis type
    if xaxis == 'lon':
        x_label = 'longitude (degrees)'
        xaxis_out = outbound_data.longitude [ ~np.isnan( outbound_data.longitude)] # np.nan_to_num() as an alternative?
        xaxis_in =  inbound_data.longitude[ ~np.isnan( inbound_data.longitude)]

    elif xaxis == 'lat':
        x_label = 'latitude (degrees)'
        xaxis_out = outbound_data.latitude[ ~np.isnan( outbound_data.latitude)]
        # min_lat = np.min( inbound_data.latitude)
        # max_lat = np.max( inbound_data.latitude)
        xaxis_in = inbound_data.latitude[ ~np.isnan( inbound_data.latitude)]
        # xaxis_in = np.linspace( max_lat, min_lat, len( xaxis_in))
        # plt.gca().invert_xaxis()

    elif xaxis == 'dist':
        xaxis_in = inbound_data.radius
        xaxis_out = - outbound_data.radius
        x_label = 'distance (km)'

    else:
        print("Error: Please Choose 'lat', 'lon', 'time', or 'distance' for the x axis")
        return

    # make plot
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()

    # plot outbound data
    # get rid of nans and resize array to get rid of overlapping data
    reflectivity = outbound_data.REFLECTIVITY.isel(time=0).isel(heading=0).transpose()

    reflectivity = reflectivity[:, range( len( xaxis_out) )]

    # print( 'number of outbound points: ' + str( len( xaxis_out)) )

    if len( xaxis_out) != 0:
        plt.pcolormesh( xaxis_out, outbound_data.height, reflectivity, cmap = color_map, vmin= -10, vmax = 50 )

    # Plot inbound data
    reflectivity = inbound_data.REFLECTIVITY.isel(time=0).isel(heading=0).transpose()

    reflectivity = reflectivity[:, range( len( xaxis_in) )]

    # print( 'number of inbound points: ' + str( len( xaxis_in)) )

    if len( xaxis_in) != 0:
        plt.pcolormesh( xaxis_in, inbound_data.height, reflectivity, cmap = color_map, vmin = -10, vmax = 50 )

    # making things prettier
    if len( xaxis_in) != 0 or len( xaxis_out) != 0:
        plt.colorbar( label="Reflectivity (dBZ)")
    plt.ylabel( 'Height from Surface (km)')
    plt.xlabel( x_label)
    plt.grid( 'on')
    # plt.gca().invert_xaxis()

    warnings.filterwarnings("default")




def plot_tdr_radial_vel( tdr_path, inbound_name, outbound_name, xaxis):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( tdr_path)
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # choose x axis type
    if xaxis == 'lon':
        x_label = 'longitude (degrees)'
        xaxis_out = outbound_data.longitude [ ~np.isnan( outbound_data.longitude)] # np.nan_to_num() as an alternative?
        xaxis_in =  inbound_data.longitude[ ~np.isnan( inbound_data.longitude)]
    elif xaxis == 'lat':
        x_label = 'latitude (degrees)'
        xaxis_out = outbound_data.latitude[ ~np.isnan( outbound_data.latitude)]
        # min_lat = np.min( inbound_data.latitude)
        # max_lat = np.max( inbound_data.latitude)
        xaxis_in = inbound_data.latitude[ ~np.isnan( inbound_data.latitude)]
        # xaxis_in = np.linspace( max_lat, min_lat, len( xaxis_in))
        # plt.gca().invert_xaxis()
    elif xaxis == 'dist':
        xaxis_in = inbound_data.radius
        xaxis_out = - outbound_data.radius
        x_label = 'distance (km)'
    else:
        print("Error: Please Choose 'lat', 'lon', 'time', or 'distance' for the x axis")
        return

    # make plot
    color_map = plt.cm.get_cmap( "RdBu").reversed()
    # set 0 as the central point of the pcolormesh!
    divnorm = colors.TwoSlopeNorm(vmin=-20, vcenter=0, vmax=20)

    # plot outbound data
    # get rid of nans and resize array to get rid of overlapping data
    vel = outbound_data.Radial_wind.isel(time=0).isel(heading=0).transpose()
    vel = vel[:, range( len( xaxis_out) )]
    plt.pcolormesh( xaxis_out, outbound_data.height, vel, cmap = color_map, norm=divnorm )

    # Plot inbound data
    vel = inbound_data.Radial_wind.isel(time=0).isel(heading=0).transpose()
    vel = vel[:, range( len( xaxis_in) )]

    # get rid of negative sign?!?
    plt.pcolormesh( xaxis_in, inbound_data.height, - vel, cmap = color_map, norm=divnorm )

    # making things prettier
    plt.colorbar( label="Radial Velocity (m/s)")
    plt.ylabel( 'Height from Surface (km)')
    plt.xlabel( x_label)
    plt.grid( 'on')
    # plt.gca().invert_xaxis()

    warnings.filterwarnings("default")



def plot_tdr_vertical_vel( tdr_path, inbound_name, outbound_name, xaxis):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( tdr_path)
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # choose x axis type
    if xaxis == 'lon':
        x_label = 'longitude (degrees)'
        xaxis_out = outbound_data.longitude [ ~np.isnan( outbound_data.longitude)] # np.nan_to_num() as an alternative?
        xaxis_in =  inbound_data.longitude[ ~np.isnan( inbound_data.longitude)]
    elif xaxis == 'lat':
        x_label = 'latitude (degrees)'
        xaxis_out = outbound_data.latitude[ ~np.isnan( outbound_data.latitude)]
        # min_lat = np.min( inbound_data.latitude)
        # max_lat = np.max( inbound_data.latitude)
        xaxis_in = inbound_data.latitude[ ~np.isnan( inbound_data.latitude)]
        # xaxis_in = np.linspace( max_lat, min_lat, len( xaxis_in))
        # plt.gca().invert_xaxis()
    elif xaxis == 'dist':
        xaxis_in = inbound_data.radius
        xaxis_out = outbound_data.radius
        x_label = 'distance (km)'
    else:
        print("Error: Please Choose 'lat', 'lon', 'time', or 'distance' for the x axis")
        return

    # make plot
    color_map = plt.cm.get_cmap( "RdBu").reversed()
    # set 0 as the central point of the pcolormesh!
    divnorm = colors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)

    # plot outbound data
    # get rid of nans and resize array to get rid of overlapping data
    vel = outbound_data.Vertical_wind.isel(time=0).isel(heading=0).transpose()
    vel = vel[:, range( len( xaxis_out) )]
    plt.pcolormesh( xaxis_out, outbound_data.height, vel, cmap = color_map, norm=divnorm )

    # Plot inbound data
    vel = inbound_data.Vertical_wind.isel(time=0).isel(heading=0).transpose()
    vel = vel[:, range( len( xaxis_in) )]

    # get rid of negative sign?!?
    plt.pcolormesh( xaxis_in, inbound_data.height, - vel, cmap = color_map, norm=divnorm )

    # making things prettier
    plt.colorbar( label="Vertical Velocity (m/s)")
    plt.ylabel( 'Height from Surface (km)')
    plt.xlabel( x_label)
    plt.grid( 'on')
    # plt.gca().invert_xaxis()

    warnings.filterwarnings("default")



def plot_all(data_path, data_file, title, index1, index2, xaxis):

    warnings.filterwarnings("ignore")

    fig = plt.figure( figsize=(20, 21))
    plt.subplot(611)
    plot_T(data_path, data_file, index1, index2, xaxis)
    plt.title( title)

    plt.subplot(612)
    plot_wvmr(data_path, data_file, index1, index2, xaxis)

    plt.subplot(613)
    plot_lsr(data_path, data_file, index1, index2, xaxis)

    plt.subplot(614)
    plot_power_ch1(data_path, data_file, index1, index2, xaxis)

    plt.subplot(615)
    plot_rh(data_path, data_file, index1, index2, xaxis)

    plt.subplot(616)
    plot_T_anomaly(data_path, data_file, index1, index2, xaxis)

    # choose x axis type
    if xaxis == 'lon':
        plt.xlabel( 'longitude (degrees)')
    elif xaxis == 'lat':
        plt.xlabel('latitude (degrees)' )
    elif xaxis == 'time':
        plt.xlabel( 'Time (UTC)' )
    elif xaxis == 'dist':
        plt.xlabel( 'Distance (km)')

    warnings.filterwarnings("default")



def plot_temps_wv(data_path, data_file, title, index1, index2, xaxis):

    warnings.filterwarnings("ignore")

    fig = plt.figure( figsize=(20, 16))
    plt.subplot(411)
    plot_T(data_path, data_file, index1, index2, xaxis)
    plt.title( title)

    plt.subplot(412)
    plot_T_anomaly(data_path, data_file, index1, index2, xaxis)

    plt.subplot(413)
    plot_wvmr(data_path, data_file, index1, index2, xaxis)

    plt.subplot(414)
    plot_rh(data_path, data_file, index1, index2, xaxis)

    # choose x axis type
    if xaxis == 'lon':
        plt.xlabel( 'longitude (degrees)')
    elif xaxis == 'lat':
        plt.xlabel('latitude (degrees)' )
    elif xaxis == 'time':
        plt.xlabel( 'Time (UTC)' )
    elif xaxis == 'dist':
        plt.xlabel( 'Distance (km)')

    warnings.filterwarnings("default")


def load_crl( path, print_files=True):
    return helper_fns.display_data_files( path, 'crl', print_files)

def load_tdr( path, print_files=True):
    return helper_fns.display_data_files( path, 'tdr', print_files)

def load_sondes( path, print_files=True):
    return helper_fns.display_data_files( path, 'dropsonde', print_files)

def load_flight_level( path, print_files=True):
    return helper_fns.display_data_files( path, 'in-situ', print_files)
