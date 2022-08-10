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
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts")
import load_in_situ_data


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
            title = "CRL and In Situ Data, TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset]

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

    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    keyList = [ '']

    for key in keyList:
        if len( cutoff_indices) == 2:
            clip_data.helper( crl_xr_in_situ[ key].values, 2)

    if len(cutoff_indices) == 2:
        ws = clip_data.helper( xr_in_situ, 'WS.d', 2)
        wd = clip_data.helper( xr_in_situ, "WD.d", 2)
        uwz = clip_data.helper( xr_in_situ, "UWZ.d", 2)
        rr = clip_data.helper( xr_in_situ, "ASfmrRainRate.1", 2)
        lat = clip_data.helper( xr_in_situ, "LATref", 2)
        lon = clip_data.helper( xr_in_situ, "LONref", 2)
        tas = clip_data.helper( xr_in_situ, "TAS.d", 2) # true air speed (m/s)
        tas2 = clip_data.helper( xr_in_situ, "TASkt.d".values, 2) # tas from knots to m/s
        height = clip_data.helper( xr_in_situ["HT.d"].values, 2 ) # height above ocean surface
        height2 = clip_data.helper( xr_in_situ[ "GPS_GeoidHt.3"].values, 2 ] # another height representation
        temp = [ float( line) for line in xr_in_situ["TTM.1"].values ] # total temperature
        temp2 = [ float( line) for line in xr_in_situ["TTMref"].values ] # raw total temperature as a reference
    elif len( cutoff_indices) == 4:
        clip_data.helper( 4)
    else:
        print( "There must only be 2 or 4 indices for trimming CRL data!")

    '''
    ws = [ float( line) for line in xr_in_situ["WS.d"].values] # convert strings to floats
    wd = [ float( line) for line in xr_in_situ["WD.d"].values]
    uwz = [ float( line) for line in xr_in_situ["UWZ.d"].values]
    rr = [ float( line) for line in xr_in_situ["ASfmrRainRate.1"].values]
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]
    tas = [float( line) for line in xr_in_situ["TAS.d"].values ] # true air speed (m/s)
    tas2 = [float( line) / 1.9438 for line in xr_in_situ["TASkt.d"].values ] # tas from knots to m/s
    height = [ float( line) for line in xr_in_situ["HT.d"].values ] # height above ocean surface
    height2 = [ float( line) for line in xr_in_situ["GPS_GeoidHt.3"].values ] # another height representation
    temp = [ float( line) for line in xr_in_situ["TTM.1"].values ] # total temperature
    temp2 = [ float( line) for line in xr_in_situ["TTMref"].values ] # raw total temperature as a reference


    # deal with annoying indices :(
    if len( cutoff_indices) == 2:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
        time1 = crl_data.time[ i1]
        time2 = crl_data.time[ i2]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()
        # use the indices for nearest time values to trim down lat and lon values
        # this prevents data overlap and makes plotting faster!
        lon = lon[ idx1.values : idx2.values]
        lat = lat[ idx1.values : idx2.values]
        ws = ws[ idx1.values : idx2.values]
        uwz = uwz[ idx1.values : idx2.values]
        rr = rr[ idx1.values : idx2.values]
        tas = tas[ idx1.values : idx2.values]
        tas2 = tas2[ idx1.values : idx2.values]
        height = height[ idx1.values : idx2.values]
        height2 = height2[ idx1.values : idx2.values]
        temp = temp[ idx1.values : idx2.values]
        temp2 = temp2[ idx1.values : idx2.values]

    elif len( cutoff_indices) == 4:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
        i3 = cutoff_indices[2]
        i4 = cutoff_indices[3]
        time1 = crl_data.time[ i1]
        time2 = crl_data.time[ i2]
        time3 = crl_data.time[ i3]
        time4 = crl_data.time[ i4]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()
        idx3 = (np.abs(float_time - time3)).argmin()
        idx4 = (np.abs(float_time - time4)).argmin()

        # use the indices for nearest time values to trim down lat and lon values
        # this prevents data overlap and makes plotting faster!
        lon = lon[ idx1.values : idx2.values] + lon[ idx3.values : idx4.values]
        lat = lat[ idx1.values : idx2.values] + lat[ idx3.values : idx4.values]
        ws = ws[ idx1.values : idx2.values] + ws[ idx3.values : idx4.values]
        uwz = uwz[ idx1.values : idx2.values] + uwz[ idx3.values : idx4.values]
        rr = rr[ idx1.values : idx2.values] + rr[ idx3.values : idx4.values]
        tas = tas[ idx1.values : idx2.values] + tas[ idx3.values : idx4.values]
        tas2 = tas2[ idx1.values : idx2.values] + tas2[ idx3.values : idx4.values]
        height = height[ idx1.values : idx2.values] + height[ idx3.values : idx4.values]
        height2 = height2[ idx1.values : idx2.values] + height2[ idx3.values : idx4.values]
        temp = temp[ idx1.values : idx2.values] + temp[ idx3.values : idx4.values]
        temp2 = temp2[ idx1.values : idx2.values] + temp2[ idx3.values : idx4.values]
    else:
        print( 'Error in number of indices! Update them in tc_metadata.py')
    '''

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

    fig = plt.figure( figsize=(14, 16) )
    ax1 = fig.add_subplot(311)

    ax1.plot( xaxis_data, ws, color='r')
    # calling this twice to set ylabel size correctly! doesn't work the first time
    ax1.set_ylabel('Tangential Wind Speed (m/s)', color='r')
    ax1.set_ylabel('Tangential Wind Speed (m/s)', color='r')
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
