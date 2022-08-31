import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
import xarray as xr
# from scipy.signal import find_peaks
from matplotlib import cm
from matplotlib.colors import ListedColormap


os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import tc_metadata
import helper_fns

# this is an attempt to get clip_old_data to load...
import sys
sys.path.append("/Users/etmu9498/research/code/scripts/in-situ-scripts/")
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts/")
import load_in_situ_data
import clip_old_data


def plot_multi_passes( tc='all'):

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    for tcname in tcname_list:
        warnings.filterwarnings("ignore")

        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return

        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')

        for dataset in range( len( metadata[ 'dates'] )):

            if metadata[ 'xlims'] [ dataset][0] > 0.0:
                axis = 'lat'
            else:
                axis = 'lon'

            # load data
            os.chdir( metadata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( metadata['dates'][dataset], metadata['crl_list'] )

            os.chdir( metadata['in_situ_path'] )
            in_situ_data = tc_metadata.choose_in_situ_date( metadata['dates'][dataset], metadata['in_situ_list'] )

            # xlims = [ metadata['xlims'][dataset][0], metadata['xlims'][dataset][1] ]
            xname = metadata['xtype'][dataset]
            xlim = metadata['xlims'][dataset]
            title = "In Situ Data, TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset]

            helper_fns.change_font_sizes()
            plot_one_cross_section( metadata['crl_path'], crl_data, metadata['in_situ_path'], in_situ_data, metadata['crl_range'][dataset], axis, xlim, xname, title)

            # save the figure!
            os.chdir( "/Users/etmu9498/research/figures/in-situ-only/")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png" )
            print( "Image " + str( dataset + 1) + " complete" )



def plot_one_cross_section( crl_path, crl_name, flight_data_path, flight_name, cutoff_indices, xaxis, xlim='none', xname='none', title='none'):

    warnings.filterwarnings("ignore")

    # load and process the data
    xr_in_situ = load_in_situ_data.load_in_situ( flight_data_path, flight_name, sample_step_size=1)

    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time

    keyList = [ 'WS.d', 'WD.d', 'UWZ.d', 'ASfmrRainRate.1', 'LATref', 'LONref',
            'TAS.d', 'HT.d', 'TTMref']
    # make an empty dict that will be filled soon!
    datatrim = {key: None for key in keyList}

    for key in keyList:
        datatrim[ key] = clip_old_data.in_situ_helper( crl_path, crl_name, cutoff_indices, xr_in_situ[ key].values, float_time)

    ws, wd, uwz = datatrim['WS.d'], datatrim['WD.d'], datatrim['UWZ.d']
    rr, lat, lon = datatrim['ASfmrRainRate.1'], datatrim['LATref'], datatrim['LONref']
    tas, height, temp2 = datatrim['TAS.d'], datatrim['HT.d'], datatrim['TTMref']

    if xaxis == "time":
        xaxis_data = float_time
        xlabel = 'Time (UTC, Hours)'
    elif xaxis == "lon":
        xaxis_data = lon
        xlabel = 'Longitude (Degrees)'
    elif xaxis == "lat":
        xaxis_data = lat
        xlabel = 'Latitude (Degrees)'
    elif xaxis == 'dist':
        print( 'add distance option to in_situ_multi panels if statement!')
        return
    else:
        print( "Please choose 'lon', 'lat', 'time', or 'dist' as a valid xaxis!")


    # find tangential wind speed and p-3 height minima! plot them to compare
    wsmin = np.nanmin( ws)
    wsmin_ind = np.nanargmin( ws)

    heightmin = np.nanmin( height)
    heightmin_ind = np.nanargmin( height)

    fig = plt.figure( figsize=(14, 16) )
    ax1 = fig.add_subplot(311)

    ax1.plot( xaxis_data, ws, color='r')
    # calling this twice to set ylabel size correctly! doesn't work the first time
    ax1.set_ylabel('Tangential Wind Speed (m/s)', color='r')

    # plot vertical lines representing wind speed and height minima
    # ax1.axvline( x= xaxis_data[ wsmin_ind], c='r')
    # ax1.axvline( x=xaxis_data[ heightmin_ind], c='y')

    ax1.xaxis.grid( )
    ax1.yaxis.grid( )
    if title != 'none':
        ax1.set_title( title)
    if xlim != 'none':
        ax1.set_xlim( xlim)

    ax2 = ax1.twinx()
    ax2.plot( xaxis_data, uwz, color='k')
    ax2.set_ylabel('Vertical Wind Speed (m/s)', color='k')
    # ax2.set_ylim([-.5, 15])

    ax3 = fig.add_subplot(312)
    ax3.plot( xaxis_data, tas, c='darkgreen')
    ax3.set_ylabel('P-3 True Air Speed (m/s)', color='g')
    ax3.xaxis.grid( )
    ax3.yaxis.grid( )
    if xlim != 'none':
        ax3.set_xlim( xlim)

    ax4 = ax3.twinx()
    ax4.plot( xaxis_data, rr, color='b')
    ax4.set_ylabel( ' SFMR Rain Rate (mm/hr)', color='b')

    ax5 = fig.add_subplot(313)
    ax5.plot( xaxis_data, height, c='y')
    ax5.set_ylabel('P-3 Height (m)', c='y')
    ax5.xaxis.grid( )
    ax5.yaxis.grid( )
    # ax5.set_ylim( [ 2500, 3300])

    # plot vertical lines representing wind speed and height minima
    # ax5.axvline( x= xaxis_data[ wsmin_ind], c='r')
    # ax5.axvline( x=xaxis_data[ heightmin_ind], c='y')


    ax6 = ax5.twinx()
    ax6.plot( xaxis_data, temp2, c='c')
    ax6.set_ylabel('Flight Level Temperature (C)', c='c')

    if xlim != 'none':
        ax5.set_xlim( xlim)

    # plt.tight_layout()
    if xname == 'lon':
        ax5.set_xlabel( "Longitude (Degrees)")
    elif xname == 'time':
        ax5.set_xlabel( "Time (UTC, Hours)")
    elif xname == 'lat':
        ax5.set_xlabel( "Latitude (Degrees)")

    warnings.filterwarnings("ignore")



def make_one_subplot( crl_path, crl_name, flight_data_path, flight_name, cutoff_indices):

    warnings.filterwarnings("ignore")

    # load and process the data

    os.chdir( flight_data_path)
    xr_in_situ = xr.open_dataset( flight_name)

    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time

    keyList = [ 'distance', 'WS.d', 'HT.d']
    # make an empty dict that will be filled soon!
    datatrim = {key: None for key in keyList}


    for key in keyList:
        datatrim[ key] = clip_old_data.in_situ_helper( crl_path, crl_name, cutoff_indices, xr_in_situ[ key].values, float_time)

    xaxis_data, ws, height = datatrim['distance'], datatrim['WS.d'], datatrim['HT.d']


    fig = plt.gcf()
    ax1 = fig.add_subplot(414)

    ax1.plot( xaxis_data, height, c='y')
    ax1.set_ylabel('P-3 Height (m)', c='y')
    # ax1.set_ylim( [ 2500, 3300])
    ax1.xaxis.grid( )
    ax1.yaxis.grid( )

    ax2 = ax1.twinx()
    ax2.plot( xaxis_data, ws, c='c')
    ax2.set_ylabel( 'Tangential Wind Speed (m/s)', c='c')

    warnings.filterwarnings("ignore")
