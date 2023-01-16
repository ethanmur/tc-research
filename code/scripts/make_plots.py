## Last edited 7/27/22


import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import xarray as xr
import warnings
import datetime
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator

os.chdir(  "/Users/etmu9498/research/code/scripts")
import helper_fns
import tc_metadata


# add a scale to try to better collocate data!
scale = 1 # 1.002
str_ticks = 10

# a helper script used in most plotting functions to determine the x axis scale variable

def x_lims_helper( axis):

    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    if axis == 'lon-str':
        xaxis = crl_data.Lon[index1:index2]
        # np.array2string( xaxis, threshold=index2, precision=4, separator=', ')
        # [ str( value) for value in xaxis]
        # xaxis.astype('|S5')
        xaxis_ticks = [ np.array2string( value) for value in xaxis]
    elif axis == 'lat-str':
        xaxis = crl_data.Lon[index1:index2]
        # np.array2string( xaxis, threshold=index2, precision=4, separator=', ')
        # [ str( value) for value in xaxis]
        # xaxis.astype('|S5')
        xaxis_ticks = [ np.array2string( value) for value in xaxis]


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

    # print the longitude axis values as strings, not floats, so data doesn't wrap on itself
    elif xaxis == 'lon-str':

        xlims = None # [ xaxis[index1], xaxis[index2 - 1 ] ]
        xaxis = range( index2 - 1)
        x_label = 'longitude (degrees)'
        return xaxis, x_label, xlims

    # print the latitude axis values as strings, not floats, so data doesn't wrap on itself
    elif xaxis == 'lat-str':
        xaxis = crl_data.Lat[index1:index2]
        xaxis = ["%.3f" % number for number in xaxis]
        xlims = None # [ xaxis[index1], xaxis[index2 - 1] ]
        x_label = 'latitude (degrees)'
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


def plot_T(data_path, data_file, index1, index2, xaxis_name, show_colorbar=True):

    warnings.filterwarnings("ignore")

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()


    # choose x axis with helper script
    xaxis, x_label, xlims = x_axis_helper( data_path, data_file, index1, index2, xaxis_name)

    # calculate temperature
    temp = crl_data.T[index1:index2, :].where( crl_data.T[index1:index2, :].values < 50).transpose()

    # plot and make things pretty
    # code to correct negative heights on some new 2022 crl data
    if np.mean( crl_data.H) > 0.0:
        plt.pcolormesh( xaxis, crl_data.H, temp, cmap = color_map, vmin=5, vmax=35 )
    else:
        plt.pcolormesh( xaxis, - crl_data.H, temp, cmap = color_map, vmin=5, vmax=35 )

    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    ax = plt.gca()
    ax.set_facecolor('k')
    if show_colorbar:
        plt.colorbar(label="Temperature ( C)")

    if xaxis_name == 'lon-str' or xaxis_name == 'lat-str':
        ax.xaxis.set_major_locator(MaxNLocator( str_ticks))

    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.grid( True)
    warnings.filterwarnings("default")


def plot_new_T( data_path, data_file, data_source = 'none', xlims=None, show_colorbar=True):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( data_path)
    new_crl = xr.open_dataset( data_file)
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()

    # print( len( new_crl.in_situ_distance))
    # print( type( new_crl.psurf_distance))

    # print( new_crl.psurf_distance)

    # print( len( new_crl.psurf_distance))

    # if np.isnan( new_crl.psurf_distance):
    #     print( 'nan')
    # else:
    #     print( 'no nan')


    if data_source == 'tdr':
        temp = new_crl.T[ 0 : len( new_crl.tdr_distance), :]
        temp = temp.where( temp < 50).transpose()
        plt.pcolormesh( new_crl.tdr_distance, - new_crl.H, temp, cmap = color_map, vmin=5, vmax=35 )
    elif data_source == 'in-situ':
        temp = new_crl.T[ 0 : len( new_crl.in_situ_distance), :]
        temp = temp.where( temp < 50).transpose()
        plt.pcolormesh( new_crl.in_situ_distance, - new_crl.H, temp, cmap = color_map, vmin=5, vmax=35 )
    elif data_source == 'in-situ-psurf':
        temp = new_crl.T[ 0 : len( new_crl.psurf_distance), :]
        temp = temp.where( temp < 50).transpose()
        plt.pcolormesh( new_crl.psurf_distance, - new_crl.H, temp, cmap = color_map, vmin=5, vmax=35 )
    else:
        print( "Please choose either 'tdr', 'in-situ', or 'in-situ-psurf' for the data_source input for the plot_new_power_ch1() function")


    if show_colorbar:
        plt.colorbar(label="Temperature ( C)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')


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
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if xaxis_name == 'lon-str' or xaxis_name == 'lat-str':
        ax.xaxis.set_major_locator(MaxNLocator( str_ticks))



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
    # code to correct negative heights on some new 2022 crl data
    if np.mean( crl_data.H) > 0.0:
        plt.pcolormesh( xaxis, crl_data.H, crl_wvmr, vmin = 0, vmax =20)
    else:
        plt.pcolormesh( xaxis, - crl_data.H, crl_wvmr, vmin = 0, vmax =20)

    plt.colorbar(label="WVMR ( g/kg)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if xaxis_name == 'lon-str' or xaxis_name == 'lat-str':
        ax.xaxis.set_major_locator(MaxNLocator( str_ticks))



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
    # code to correct negative heights on some new 2022 crl data
    if np.mean( crl_data.H) > 0.0:
        plt.pcolormesh( xaxis, crl_data.H, crl_lsr)
    else:
        plt.pcolormesh( xaxis, - crl_data.H, crl_lsr)

    plt.colorbar(label="LSR ( unitless)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if xaxis_name == 'lon-str' or xaxis_name == 'lat-str':
        ax.xaxis.set_major_locator(MaxNLocator( str_ticks))



def plot_power_ch1(data_path, data_file, index1, index2, xaxis_name, cutoff=-30, show_colorbar=True):

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

    # plot things

    # looking at different colormaps
    # color_map = nclcmaps.cmap( 'MPL_brg'  "MPL_RdGy" 'ncview_default')
    # color_map = plt.cm.get_cmap( "RdYlBu").reversed()
    # color_map = 'viridis'
    # color_map = plt.cm.get_cmap( "Spectral").reversed()
    # code to correct negative heights on some new 2022 crl data
    if np.mean( crl_data.H) > 0.0:
        plt.pcolormesh(  xaxis, crl_data.H, crl_pch1, vmin = cutoff, vmax =-10)
    else:
        plt.pcolormesh(  xaxis, - crl_data.H, crl_pch1, vmin = cutoff, vmax =-10)

    if show_colorbar:
        plt.colorbar(label="CRL Return Power (dBz)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if xaxis_name == 'lon-str' or xaxis_name == 'lat-str':
        ax.xaxis.set_major_locator(MaxNLocator( str_ticks))

    warnings.filterwarnings("default")


def plot_new_power_ch1( data_path, data_file, data_source='none', cutoff=-30, xlims=None, show_colorbar=True):
    warnings.filterwarnings("ignore")


    # get data
    os.chdir( data_path)
    new_crl = xr.open_dataset( data_file)
    # manipulate data
    step1 = 10 * np.log10( new_crl.P_ch1 )
    step2 = step1.where( step1.values > cutoff)

    if data_source == 'tdr':
        # cut off the last ind to fit the x axis!
        crl_pch1 = step2[0 : len( new_crl.tdr_distance) , : ].transpose()
        plt.pcolormesh(  new_crl.tdr_distance, - new_crl.H, crl_pch1, vmin = cutoff, vmax =-10)
    elif data_source == 'in-situ':
        # cut off the last ind to fit the x axis!
        crl_pch1 = step2[0 : len( new_crl.in_situ_distance) , : ].transpose()

        # check if there's a nan in the distance array
        if np.isnan( new_crl.in_situ_distance).any() :
            print( 'nan!')

        plt.pcolormesh(  new_crl.in_situ_distance, - new_crl.H, crl_pch1, vmin = cutoff, vmax =-10)
    else:
        print( "Please choose either 'tdr' or 'in-situ' for the data_source input for the plot_new_power_ch1() function")

    if show_colorbar:
        plt.colorbar(label="Backscattered Ch 1 power ( dBz)")
    plt.ylabel( 'Height (km)')
    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')



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

    '''
    print( 'temp: ' + str( type( temp)))
    print( 'temp: ' + str( np.shape( temp)))
    print( 'wvmr: ' + str( type( wvmr)))
    print( 'wvmr: ' + str( np.shape( wvmr)))
    print( 'pressure: ' + str( type( pressure)))
    print( 'pressure: ' + str( np.shape( pressure)))
    print( 'es: ' + str( type( e_s)))
    print( 'es: ' + str( np.shape( e_s)))
    print( 'rh: ' + str( type( crl_rh)))
    print( 'rh: ' + str( np.shape( crl_rh)))
    '''

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H, crl_rh)
    plt.colorbar(label="RH ( %)")
    plt.ylabel( 'Height (km)')

    # unique fix for rh code axes, they need to be flipped for some reason to get plots to match :/
    if xaxis_name == 'lon' or xaxis_name == 'lat':
        plt.gca().invert_xaxis()

    # plt.xlabel( x_label)
    if xlims:
        plt.xlim( xlims )
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if xaxis_name == 'lon-str' or xaxis_name == 'lat-str':
        ax.xaxis.set_major_locator(MaxNLocator( str_ticks))



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




def plot_tdr( tdr_path, inbound_name, outbound_name, xaxis):
    # the function I'm using to plot tdr reflectivity

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
        xaxis_in = - inbound_data.radius
        xaxis_out = outbound_data.radius
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
    # plt.xlabel( x_label)
    plt.grid( 'on')
    # plt.gca().invert_xaxis()

    warnings.filterwarnings("default")


# This function plots the reflectivity data from the new, combined tdr files I
# created for every good eye case
def plot_new_tdr( tdr_path, tdr_name, xaxis='dist'):
    warnings.filterwarnings("ignore")

    # get data
    os.chdir( tdr_path)
    tdr_data = xr.open_dataset( tdr_name)

    # choose x axis type
    if xaxis == 'lon':
        x_label = 'longitude (degrees)'
        xaxis = tdr_data.longitude[ ~np.isnan( tdr_data.longitude)]
        '''
        # getting rid of nans on either end of latitude data to avoid plotting errors
        # not entirely sure why I'm doing it this way... wouldn't just sorting for nans work?
        # make a list of indices for latitude values
        xinds = range( len( tdr_data.longitude))
        # turn all indices with a Nan value into a 0
        # a value other than 0 can be used here too, maybe something really big
        xaxis = np.where( np.isnan( tdr_data.longitude), 0, xinds)
        # get rid of all zeros
        # this might unfairly get rid of the first data point (aka ind 0), but that's usually a nan anyways!
        xaxis_ind = xaxis[ xaxis != 0]
        # only select non nan latitude values
        xaxis = tdr_data.longitude[ xaxis_ind]
        '''

    elif xaxis == 'lat':
        x_label = 'latitude (degrees)'
        xaxis = tdr_data.latitude[ ~np.isnan( tdr_data.latitude)]

    elif xaxis == 'dist':
        xaxis = tdr_data.distance
        x_label = 'distance (km)'
    else:
        print("Error: Please Choose 'lat', 'lon', or 'dist' for the x axis")
        return

    # make plot
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()

    # plot data
    # get rid of nans and resize array to get rid of overlapping data
    # also, no need to use .transpose() because that was already done when making the datasets!
    reflectivity = tdr_data.REFLECTIVITY[ :, 0:len( xaxis)]
    # refl = tdr_data.REFLECTIVITY[:, lat_no_nan_ind] # another way?
    plt.pcolormesh( xaxis, tdr_data.height, reflectivity, cmap = color_map, vmin = -10, vmax = 50 )

    # making things prettier
    if len( xaxis) != 0 or len( xaxis) != 0:
        plt.colorbar( label="Reflectivity (dBZ)")
    plt.ylabel( 'Height (Km)')
    plt.grid( 'on')
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


def plot_new_radial_vel( tdr_path, tdr_name, xaxis='dist'):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( tdr_path)
    tdr_data = xr.open_dataset( tdr_name)

    # choose x axis type
    if xaxis == 'lon':
        x_label = 'longitude (degrees)'
        xaxis = tdr_data.longitude[ ~np.isnan( tdr_data.longitude)]
    elif xaxis == 'lat':
        x_label = 'latitude (degrees)'
        xaxis = tdr_data.latitude[ ~np.isnan( tdr_data.latitude)]

    elif xaxis == 'dist':
        xaxis = tdr_data.distance
        x_label = 'distance (km)'
    else:
        print("Error: Please Choose 'lat', 'lon', or 'dist' for the x axis")
        return

    # make plot
    color_map = plt.cm.get_cmap( "RdBu").reversed()
    # set 0 as the central point of the pcolormesh!
    divnorm = colors.TwoSlopeNorm(vmin=-20, vcenter=0, vmax=20)

    # plot all data
    # get rid of nans and resize array to get rid of overlapping data
    vel = tdr_data.Radial_wind[ :, 0:len( xaxis)] # .transpose()
    plt.pcolormesh( xaxis, tdr_data.height, vel, cmap = color_map, norm=divnorm )

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






def plot_new_vertical_vel( tdr_path, tdr_name, xaxis='dist'):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( tdr_path)
    tdr_data = xr.open_dataset( tdr_name)

    # choose x axis type
    if xaxis == 'lon':
        x_label = 'longitude (degrees)'
        xaxis = tdr_data.longitude[ ~np.isnan( tdr_data.longitude)]
    elif xaxis == 'lat':
        x_label = 'latitude (degrees)'
        xaxis = tdr_data.latitude[ ~np.isnan( tdr_data.latitude)]

    elif xaxis == 'dist':
        xaxis = tdr_data.distance
        x_label = 'distance (km)'
    else:
        print("Error: Please Choose 'lat', 'lon', or 'dist' for the x axis")
        return

    # make plot
    color_map = plt.cm.get_cmap( "RdBu").reversed()
    # set 0 as the central point of the pcolormesh!
    divnorm = colors.TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)

    # plot all data
    # get rid of nans and resize array to get rid of overlapping data
    vel = tdr_data.Vertical_wind[ :, 0:len( xaxis)] # .transpose()
    plt.pcolormesh( xaxis, tdr_data.height, vel, cmap = color_map, norm=divnorm )

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

def plot_some(data_path, data_file, title, index1, index2, xaxis):

    warnings.filterwarnings("ignore")

    fig = plt.figure( figsize=(20, 8))
    plt.subplot(211)
    plot_T(data_path, data_file, index1, index2, xaxis)
    plt.title( title)

    plt.subplot(212)
    plot_power_ch1(data_path, data_file, index1, index2, xaxis)

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


def plot_full_dataset_one_day( data_path, data_file, xaxis, title=None):
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    index1 = 0
    index2 = len( crl_data.time) - 1

    warnings.filterwarnings("ignore")
    fig = plt.figure( figsize=(22, 24))

    plt.subplot(611)
    plot_T(data_path, data_file, index1, index2, xaxis)
    if title:
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
        plt.xlabel( 'Time (UTC, Hours)' )
    elif xaxis == 'dist':
        plt.xlabel( 'Distance (km)')



def plot_full_datasets( tc='all', xaxis='time'):
    if tc == 'all':
        tcname_list = ['fred', 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]
    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        tcdata = tc_metadata.plot_all_crl_data( tcname)

        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['dates'] ))+ '\n')
        # load data

        # plot each day's dataset for the tc
        for counter in range( len( tcdata[ 'dates'])):
            # get the correct name of this day's dataset
            crl_path = tcdata[ 'crl_path']
            crl_name = tc_metadata.choose_crl_date( tcdata[ 'dates'][counter], tcdata[ 'crl_list'])

            title = "Full CRL Dataset, TC " + tcdata['tc_name'] + ", " + tcdata[ 'dates'][counter]
            plot_full_dataset_one_day( crl_path, crl_name, xaxis, title=title)

            os.chdir( "/Users/etmu9498/research/figures/all-data/")
            plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( counter + 1) + " saved\n" )



def load_crl( path, print_files=True):
    return helper_fns.display_data_files( path, 'crl', print_files)

def load_tdr( path, print_files=True):
    return helper_fns.display_data_files( path, 'tdr', print_files)

def load_sondes( path, print_files=True):
    return helper_fns.display_data_files( path, 'dropsonde', print_files)

def load_flight_level( path, print_files=True):
    return helper_fns.display_data_files( path, 'in-situ', print_files)
