## Last edited 8/29/22
# these functions only work with the height corrected matrices found under data/crl-new-matrices

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

def x_axis_helper( data_path, data_file, xaxis):

    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    # choose x axis type
    if xaxis == 'lon':
        xaxis = crl_data.Lon
        return xaxis
    elif xaxis == 'lat':
        xaxis = crl_data.Lat
        return xaxis
    elif xaxis == 'time':
        xaxis = crl_data.time
        return xaxis
    elif xaxis == 'tdr-dist':
        xaxis = crl_data.tdr_distance
        return xaxis
    elif xaxis == 'in-situ-dist':
        xaxis = crl_data.in_situ_distance
        return xaxis
    elif xaxis == 'rmw':
        # xaxis = [ str( value) for value in crl_data.rmw]
        xaxis = crl_data.rmw
        return xaxis
    elif xaxis == 'rmw-negatives':
        # xaxis = [ str( value) for value in crl_data.rmw]

        xaxis = np.ma.masked_invalid( crl_data.rmw_negatives.values)
        return xaxis
    else:
        print("Error: Please Choose 'lat', 'lon', 'time', 'in-situ-dist' or 'tdr-dist' for the x axis")
        return None


def plot_T( data_path, data_file, xaxis = 'in-situ-dist', xlims=None, show_colorbar=True, grids=True):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( data_path)
    new_crl = xr.open_dataset( data_file)
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()

    xaxis = x_axis_helper( data_path, data_file, xaxis)
    plt.pcolormesh( xaxis, - new_crl.H_new, new_crl.T_new.transpose(), cmap = color_map, vmin=5, vmax=35 )

    if show_colorbar:
        plt.colorbar(label="Temperature ( C)")
    plt.ylabel( 'Height (km)')

    if grids:
        plt.grid( 'on')

    ax = plt.gca()
    ax.set_facecolor('k')
    warnings.filterwarnings("default")


def plot_wvmr(data_path, data_file, xaxis):
    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis = x_axis_helper( data_path, data_file, xaxis)

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H_new, crl_data.wvmr_new.transpose(), vmin = 0, vmax =20, cmap="viridis")
    plt.colorbar(label="WVMR ( g/kg)")
    plt.ylabel( 'Height (km)')
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')


def plot_lsr(data_path, data_file, xaxis):
    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)
    # choose x axis with helper script
    xaxis = x_axis_helper( data_path, data_file, xaxis)

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H_new, crl_data.lsr_new.transpose() )
    plt.colorbar(label="LSR ( unitless)")
    plt.ylabel( 'Height (km)')
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')



def plot_power_ch1(data_path, data_file, xaxis, grids=True):
    warnings.filterwarnings("ignore")
    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)


    # choose x axis with helper script
    xaxis = x_axis_helper( data_path, data_file, xaxis)

    # plot things
    plt.pcolormesh(  xaxis, - crl_data.H_new, crl_data.power_new.transpose(), vmin = -30, vmax =-10, cmap='viridis')
    plt.ylabel( 'Height (km)')
    plt.colorbar(label="Backscattered Power ( dBz)")

    if grids:
        plt.grid( 'on')

    ax = plt.gca()
    ax.set_facecolor('k')
    warnings.filterwarnings("default")


def plot_rh(data_path, data_file, xaxis):
    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis = x_axis_helper( data_path, data_file, xaxis)

    # plot things
    plt.pcolormesh( xaxis, - crl_data.H_new, crl_data.rh_new.transpose() , cmap = 'Blues', vmin = 60, vmax = 100 )
    plt.colorbar(label="RH ( %)")
    plt.ylabel( 'Height (km)')

    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')

def plot_T_anomaly(data_path, data_file, xaxis):
    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # choose x axis with helper script
    xaxis = x_axis_helper( data_path, data_file, xaxis)

    # plot things
    divnorm = colors.TwoSlopeNorm(vmin=-5, vcenter=0, vmax=15)
    plt.pcolormesh( xaxis, - crl_data.H_new, crl_data.T_anomaly.transpose() , cmap = 'bwr', norm=divnorm )
    plt.colorbar(label="Temperature Anomaly ( C)")
    plt.ylabel( 'Height (km)')

    # unique fix for rh code axes, they need to be flipped for some reason to get plots to match :/
    # if xaxis_name == 'lon' or xaxis_name == 'lat':
    #     plt.gca().invert_xaxis()

    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')


# plot regions of rainfall in green to highlight these spots! And to eventually use them
# for statistics
def plot_rainfall( data_path, data_file, xaxis):

    plt.figure( figsize=(16, 6))

    # make a base power plot
    plot_power_ch1(data_path, data_file, xaxis)
    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    # run and plot rainfall locations


def plot_all(data_path, data_file, xaxis, padding, title=None):

    # get data
    os.chdir( data_path)
    crl_data = xr.open_dataset( data_file)

    warnings.filterwarnings("ignore")

    fig = plt.figure( figsize=(20, 26))
    plt.subplot(611)
    if title:
        plt.title( title)
    plot_T(data_path, data_file, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ - padding, padding])

    plt.subplot(612)
    plot_wvmr(data_path, data_file, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ - padding, padding])

    plt.subplot(613)
    plot_lsr(data_path, data_file, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ - padding, padding])

    plt.subplot(614)
    plot_power_ch1(data_path, data_file, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ - padding, padding])

    plt.subplot(615)
    plot_rh(data_path, data_file, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ - padding, padding])

    plt.subplot(616)
    plot_T_anomaly(data_path, data_file, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ - padding, padding])

    # choose x axis type
    if xaxis == 'lon':
        plt.xlabel( 'longitude (degrees)')
    elif xaxis == 'lat':
        plt.xlabel('latitude (degrees)' )
    elif xaxis == 'time':
        plt.xlabel( 'Time (UTC)' )
    elif xaxis == 'tdr-dist' or xaxis == 'in-situ-dist':
        plt.xlabel( 'Distance (km)')

    warnings.filterwarnings("default")




def plot_some(data_path, data_file, xaxis):

    warnings.filterwarnings("ignore")

    fig = plt.figure( figsize=(20, 8))
    plt.subplot(211)
    plot_T(data_path, data_file, xaxis)

    plt.subplot(212)
    plot_power_ch1(data_path, data_file, xaxis)

    # choose x axis type
    if xaxis == 'lon':
        plt.xlabel( 'longitude (degrees)')
    elif xaxis == 'lat':
        plt.xlabel('latitude (degrees)' )
    elif xaxis == 'time':
        plt.xlabel( 'Time (UTC)' )
    elif xaxis == 'tdr-dist' or xaxis == 'in-situ-dist':
        plt.xlabel( 'Distance (km)')

    warnings.filterwarnings("default")

def plot_wv_power_t( crl_path, crl_name, xaxis, padding, title):
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    warnings.filterwarnings("ignore")

    fig = plt.figure( figsize=(20, 14))

    plt.subplot(311)
    plt.title( title)
    plot_wvmr(crl_path, crl_name, xaxis)
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ padding, - padding])

    plt.subplot(312)
    plot_power_ch1( crl_path, crl_name, xaxis)
    # plt.ylim( [ -0.1, .5])
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ padding, - padding])

    plt.subplot(313)
    plot_T(crl_path, crl_name, xaxis)
    # plt.ylim( [ -0.1, .5])
    plt.ylim( [ 0, crl_data.H_max + .1])
    plt.xlim( [ padding, - padding])
    plt.xlabel( "Distance from TC Center (km)")





def plot_full_datasets( tc='all', xaxis='in-situ-dist', padding=250):
    if tc == 'all':
        tcname_list = [ 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # run this now for nice font sizes later!
    helper_fns.change_font_sizes()

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        tcdata = tc_metadata.all_data( tcname)

        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['dates'] ))+ '\n')
        # load data

        # plot each day's dataset for the tc
        for counter in range( len( tcdata[ 'dates'])):
            # get the correct name of this day's dataset
            crl_path = "/Users/etmu9498/research/data/crl-new-matrices"
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, counter)

            title = "Height Corrected CRL Data, TC " + tcdata['tc_name'] + ", " + tcdata[ 'dates'][counter] + ' Eye Pass ' + tcdata[ 'eye_pass'][counter]

            # plot_all( crl_path, crl_name, xaxis, padding, title=title)
            plot_wv_power_t( crl_path, crl_name, xaxis, padding, title=title)

            print( "Plot " + str( counter + 1) + " created\n" )

            # os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/T-boundary")
            os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/wv-temp-power")
            '''
            if padding < 250:
                os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/zoom-in-test")
            else:
                os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/zoom-out")
            '''

            plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( counter + 1) + " saved\n" )
