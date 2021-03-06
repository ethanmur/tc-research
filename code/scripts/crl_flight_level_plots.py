import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
import xarray as xr
from scipy.signal import find_peaks
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap


os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import plot_in_situ
import tc_metadata



def flight_level_colorbar_auto( tc='all', variable_list=['rr']):

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    for tcname in tcname_list:
        warnings.filterwarnings("ignore")

        tcdata = tc_metadata.choose_data_cloud_tops_good_data( tcname)
        if tcdata == 'selected TC name is not yet implemented':
            print( tcdata)
            return

        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['xlims'] ))+ '\n')

        for counter in range( len( tcdata[ 'dates'] )):

            if tcdata[ 'xlims'] [ counter][0] > 0.0:
                axis = 'lat'
            else:
                axis = 'lon'

            # load data
            os.chdir( tcdata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

            os.chdir( tcdata['in_situ_path'] )
            in_situ_data = tc_metadata.choose_in_situ_date( tcdata['dates'][counter], tcdata['in_situ_list'] )

            # make a crl plot with variable colorbars for all tc cases
            '''
            if len ( tcdata['crl_range'][counter]) == 2:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][1]
            elif len( tcdata['crl_range'][counter]) == 4:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][3]
            else:
                print( 'Error in number of indices! Update them in tc_metadata.py')
            '''

            xlims = [ tcdata['xlims'][counter][0], tcdata['xlims'][counter][1] ]

            title = "CRL and In Situ Data, TC " + tcdata['tc_name'] + ", " + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter]
            flight_level_colorbar( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_data, tcdata['crl_range'][counter], axis, variable_list, xlims, title)

            # case to plot and save lat or lon values vs time!
            # helpful for determining if i1 and i2 are accurate
            # comment out some of the matching code above to get this to work
            # title = "In Situ Lat Lon Peaks, TC " + tcdata['tc_name'] + ", " + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter]
            # plot_lat_lon( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_data, i1, i2, axis, variable_list, xlims, title)
            # os.chdir( "/Users/etmu9498/research/figures/lat-lons/")

            # save the figure!
            os.chdir( "/Users/etmu9498/research/figures/in-situ-colorbars/")
            fig_name = ''
            for in_situ_var_name in variable_list:
                fig_name += in_situ_var_name
                fig_name += '-'

            plt.savefig( fig_name + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png", bbox_inches='tight')
            print( "Image " + str( counter + 1) + " complete" )



def flight_level_lines_auto( tc='all', variable_list=['uwz', 'ws', 'rr']):

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    for tcname in tcname_list:
        warnings.filterwarnings("ignore")

        tcdata = tc_metadata.choose_data_cloud_tops_good_data( tcname)
        if tcdata == 'selected TC name is not yet implemented':
            print( tcdata)
            return

        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['xlims'] ))+ '\n')

        for counter in range( len( tcdata[ 'dates'] )):

            if tcdata[ 'xlims'] [ counter][0] > 0.0:
                axis = 'lat'
            else:
                axis = 'lon'

            # load data
            os.chdir( tcdata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

            os.chdir( tcdata['in_situ_path'] )
            in_situ_data = tc_metadata.choose_in_situ_date( tcdata['dates'][counter], tcdata['in_situ_list'] )

            '''
            # make a crl plot with variable colorbars for all tc cases
            if len ( tcdata['crl_range'][counter]) == 2:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][1]
            elif len( tcdata['crl_range'][counter]) == 4:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][3]
            else:
                print( 'Error in number of indices! Update them in tc_metadata.py')
            '''

            xlims = [ tcdata['xlims'][counter][0], tcdata['xlims'][counter][1] ]
            title = "CRL and In Situ Data, TC " + tcdata['tc_name'] + ", " + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter]
            flight_level_lines( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_data, tcdata['crl_range'][counter], axis, variable_list, xlims, title)

            # save the figure!
            os.chdir( "/Users/etmu9498/research/figures/in-situ-lines/")
            plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
            print( "Image " + str( counter + 1) + " complete" )




def flight_level_colorbar( crl_path, crl_name, flight_data_path, flight_name, cutoff_indices, xaxis, variable_list, xlim='none', title=False):

    warnings.filterwarnings("ignore")

    # define figure
    # auto scale the size of the figure depending on the number of in situ datasets selected
    fig_len = 16 + 2 * (len( variable_list) - 1 )
    plt.figure( figsize=( fig_len, 4))
    # update font sizes for this plot!
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

    # load and process the data
    xr_in_situ = plot_in_situ.load_in_situ( flight_data_path, flight_name, sample_step_size=10)

    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time
    ws = [ float( line) for line in xr_in_situ["WS.d"].values] # convert strings to floats
    wd = [ float( line) for line in xr_in_situ["WD.d"].values]
    uwz = [ float( line) for line in xr_in_situ["UWZ.d"].values]
    rr = [ float( line) for line in xr_in_situ["ASfmrRainRate.1"].values]
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]

    # deal with annoying indices :(
    if len( cutoff_indices) == 2:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
    elif len( cutoff_indices) == 4:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[3]
    else:
        print( 'Error in number of indices! Update them in tc_metadata.py')


    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)
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

    if xaxis == "time":
        xaxis_data = float_time
        xlabel = 'Time (UTC, Hours)'
    elif xaxis == "lon":
        xaxis_data = lon
        xlabel = 'Longitude (Degrees)'
    elif xaxis == "lat":
        xaxis_data = lat
        xlabel = 'Latitude (Degrees)'
    else:
        print( "Please choose 'lon', 'lat', or 'time' as a valid xaxis!")

    # make base crl plot
    # make_plots.plot_T( crl_path, crl_name, i1, i2, xaxis)
    make_plots.plot_power_ch1( crl_path, crl_name, i1, i2, xaxis)

    crl_fig = plt.gcf()
    crl_axes = plt.gca()

    # add titles etc
    if title:
        crl_axes.set_title( title)

    crl_axes.set_xlabel( xlabel)

    if xlim != 'none':
        plt.xlim( xlim)

    in_situ_height = 3.50 # assumed height of the P3: maybe make this variable?
    offset = .1 # thickness of the flight level color mesh
    line_offset = .03 # a scalar that changes the location of the axhline
    linewidth= horiz_line_width = 2

    # add extra space on the y axis depending on the number of colorbars
    # first offset is for original colorbar, 2 * offset is for following colorbars
    plt.ylim( [0, in_situ_height + offset + 2.25 * offset * ( len( variable_list) - 1 )   ])

    # define colormaps for possible in situ datasets
    ws_cmap = 'YlOrRd' # wind speed (m/s)
    uwz_cmap = 'seismic'  # vertical velocity (m/s)
    rr_cmap = 'Blues' # rain rate from the asfmr instrument ()

    # the for loop lets us stack multiple colorbars on top of one another! and to properly scale the colorbar locatins
    for in_situ_i in range( len( variable_list)):

        y_mesh_size = [in_situ_height + 2.25 * in_situ_i * offset,
                in_situ_height + offset + 2.25 * in_situ_i * offset]

        x_cbar_0, y_cbar_0, x_cbar_size, y_cbar_size = .85 + .06 * in_situ_i, 0.125, 0.005, 0.755

        axh_height = in_situ_height - line_offset + 2 * in_situ_i * offset

        # tangential velocity case
        if variable_list[ in_situ_i] == 'ws':
            w_min = np.nanmin( ws)
            w_max = np.nanmax( ws)
            # plot in situ data above crl data in blank spot at flight level where crl data is missing
            crl_axes.pcolormesh( xaxis_data, y_mesh_size, [ws, ws], cmap = ws_cmap, vmin= w_min, vmax= w_max)
            # make a black line separating in situ and crl data
            crl_axes.axhline( y= axh_height, c='k', linewidth= horiz_line_width)

            # trying to place the colorbar closer to the figure
            # the isolated numbers determine the position of the initial colorbar, and the
            # values * counter determine the following colorbar locations
            ax_cbar = plt.gcf().add_axes([x_cbar_0, y_cbar_0, x_cbar_size, y_cbar_size])

            map = mpl.cm.ScalarMappable(cmap= ws_cmap, norm=mpl.colors.Normalize( vmin= w_min, vmax= w_max))
            plt.colorbar( mappable=map, cax=ax_cbar, orientation='vertical', label='Tangential Wind Speed (m/s)')

        # vertical velocity case
        elif variable_list[ in_situ_i] == 'uwz':
            w_min = np.nanmin( uwz)
            w_max = np.nanmax( uwz)
            # plot in situ data above crl data in blank spot at flight level where crl data is missing
            crl_axes.pcolormesh( xaxis_data, y_mesh_size, [uwz, uwz], cmap = uwz_cmap, vmin= w_min, vmax= w_max)
            # make a black line separating in situ and crl data
            crl_axes.axhline( y= axh_height, c='k', linewidth= horiz_line_width)

            # trying to place the colorbar closer to the figure
            ax_cbar = plt.gcf().add_axes([x_cbar_0, y_cbar_0, x_cbar_size, y_cbar_size])

            map = mpl.cm.ScalarMappable(cmap= uwz_cmap, norm=mpl.colors.Normalize( vmin= w_min, vmax= w_max))
            plt.colorbar( mappable=map, cax=ax_cbar, orientation='vertical', label='Vertical Wind Speed (m/s)')

        # tangential velocity case
        elif variable_list[ in_situ_i] == 'rr':
            w_min = np.nanmin( rr)
            w_max = np.nanmax( rr)
            # plot in situ data above crl data in blank spot at flight level where crl data is missing
            crl_axes.pcolormesh( xaxis_data, y_mesh_size, [rr, rr], cmap = rr_cmap, vmin= w_min, vmax= w_max)
            # make a black line separating in situ and crl data
            crl_axes.axhline( y= axh_height, c='k', linewidth= horiz_line_width)

            # trying to place the colorbar closer to the figure
            # the isolated numbers determine the position of the initial colorbar, and the
            # values * counter determine the following colorbar locations
            ax_cbar = plt.gcf().add_axes([x_cbar_0, y_cbar_0, x_cbar_size, y_cbar_size])

            map = mpl.cm.ScalarMappable(cmap= rr_cmap, norm=mpl.colors.Normalize( vmin= w_min, vmax= w_max))
            plt.colorbar( mappable=map, cax=ax_cbar, orientation='vertical', label=' SFMR Rain Rate (mm/hr)')

        else:
            print( "Please enter a valid flight level variable. These include ['ws'], ['uwz'], and ['rr'].")
            return

    warnings.filterwarnings("ignore")




def flight_level_lines( crl_path, crl_name, flight_data_path, flight_name, cutoff_indices, xaxis, xlim='none', title=False):

    warnings.filterwarnings("ignore")

    # define figure
    # auto scale the size of the figure depending on the number of in situ datasets selected
    # fig, ax = plt.subplots(2, 2, figsize=( fig_len, 8), gridspec_kw= { 'width_ratios': [1, .1]})
    plt.figure( figsize=( 16, 12))

    # update font sizes for this plot!
    TINY_SIZE = 12
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=TINY_SIZE)    # legend fontsize

    # load and process the data
    xr_in_situ = plot_in_situ.load_in_situ( flight_data_path, flight_name, sample_step_size=1) # 10

    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time
    ws = [ float( line) for line in xr_in_situ["WS.d"].values] # convert strings to floats
    wd = [ float( line) for line in xr_in_situ["WD.d"].values]
    uwz = [ float( line) for line in xr_in_situ["UWZ.d"].values]
    rr = [ float( line) for line in xr_in_situ["ASfmrRainRate.1"].values]
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]


    print( len( float_time))
    print( len( ws))
    print( len( uwz))
    print( len( rr))
    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

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
        float_time = float_time[ idx1.values : idx2.values]
        lon = lon[ idx1.values : idx2.values]
        lat = lat[ idx1.values : idx2.values]
        ws = ws[ idx1.values : idx2.values]
        uwz = uwz[ idx1.values : idx2.values]
        rr = rr[ idx1.values : idx2.values]

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

        # this code resets i2 for crl plotting below
        i2 = cutoff_indices[3]
    else:
        print( 'Error in number of indices! Update them in tc_metadata.py')


    if xaxis == "time":
        xaxis_data = float_time
        xlabel = 'Time (UTC, Hours)'
    elif xaxis == "lon":
        xaxis_data = lon
        xlabel = 'Longitude (Degrees)'
    elif xaxis == "lat":
        xaxis_data = lat
        xlabel = 'Latitude (Degrees)'
    else:
        print( "Please choose 'lon', 'lat', or 'time' as a valid xaxis!")

    plt.subplot(311)

    """
    # make sure the x axis data is a numpy array for proper plotting
    xaxis_data = np.array(xaxis_data)

    # find wind speed max points
    peaks = find_peaks( ws, prominence= 3, width=5, height= 15) # wlen=100

    # find minima to figure out where the center of the tc eye is
    # flip the wind speed dataset to use find peaks to get minima
    ws_flipped = [ -1 * val for val in ws]
    troughs = find_peaks( ws_flipped, prominence= 3, width=5, height= - 25)

    # find th, the height of the lowest trough (center of the eye)
    # th = trough height
    # - sign is to undo - sign above
    th = np.min( - troughs[1]['peak_heights'])
    # find the x position of the trough height
    # use xaxis_data[ ...] to go from an index to an x position
    # ti = trough index, tx = trough position
    ti = np.where( - troughs[1]['peak_heights'] == th)
    tx = xaxis_data[ troughs[0][ti]]
    # turn tx from a one element list into just that element
    tx = tx[0]

    # find the nearest left and right peaks to the trough tx!
    # left case:
    # find the peaks to the left of the tc eye center
    l_peaks = xaxis_data[ peaks [0]] [np.where( xaxis_data[peaks[0]] < tx)]
    # full list case (there are values to the left of the eye)
    if np.any( l_peaks):
        # lp(i, x, h) = left peak (index, position, height)
        lpi = (np.abs( l_peaks - tx )).argmin()
        lpx = l_peaks[ lpi]
        lph = ws[ np.where( xaxis_data == lpx) [0][0] ]
    # empty list case
    else:
        # just guess where the eyewall is lol... this needs to be fixed
        lpx = xaxis_data[ np.where( xaxis_data == tx)[0][0] - 50]
        lph = ws[ np.where( xaxis_data == lpx) [0][0] ]

    # right case
    r_peaks = xaxis_data[ peaks [0]] [np.where( xaxis_data[peaks[0]] >= tx)]
    # full list case
    if np.any( r_peaks):
        rpi = (np.abs( r_peaks - tx )).argmin()
        rpx = r_peaks[ rpi]
        rph = ws[ np.where( xaxis_data == rpx) [0][0] ]
    # empty list case
    else:
        rpx = xaxis_data[ np.where( xaxis_data == tx)[0][0] + 50]
        rph = ws[ np.where( xaxis_data == rpx) [0][0] ]
    # plot peaks and troughs
    # plt.scatter( xaxis_data[ peaks[0]], peaks[1]['peak_heights'], s=30, c='k', marker='x')
    # plt.scatter( xaxis_data[ troughs[0]], - troughs[1]['peak_heights'], s=30, c='b', marker='x')

    # plot the central peaks and trough in a different color and shape!
    # plt.scatter( [lpx, rpx, tx], [lph, rph, th], s=40, c='g', marker='s')


    """

    # for y_var in variable_list:
    plt.plot( xaxis_data, ws, c='r', label='Tangential Wind Speed (m/s)')
    plt.plot( xaxis_data, uwz, c='k', label='Vertical Wind Speed (m/s)')
    # plt.plot( xaxis_data, rr, c='b', label=' SFMR Rain Rate (mm/hr)')

    # add an empty colorbar to make everything fit in line... kinda a
    # messy solution but it's ok for now!
    viridis = cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    white = np.array([ 1, 1, 1, 1])
    newcolors[:, :] = white
    white_cmap = ListedColormap(newcolors)
    map = mpl.cm.ScalarMappable(cmap= white_cmap, norm=mpl.colors.Normalize( vmin= 0, vmax= 1))
    cbar = plt.colorbar( mappable= map)
    cbar.set_ticks([])
    cbar.outline.set_visible(False)


    plt.legend(loc='upper left')

    # add titles etc
    if title:
        plt.title( title)

    if xlim != 'none':
        plt.xlim( xlim)

    plt.grid('on')

    # make base crl plots
    plt.subplot(312)
    make_plots.plot_T( crl_path, crl_name, i1, i2, xaxis, show_colorbar=True)

    if xlim != 'none':
        plt.xlim( xlim)

    # plt.sca(ax[1, 0])  # set the current axes instance
    plt.subplot(313)
    make_plots.plot_power_ch1( crl_path, crl_name, i1, i2, xaxis, show_colorbar=True)

    plt.xlabel( xlabel)
    if xlim != 'none':
        plt.xlim( xlim)

    warnings.filterwarnings("ignore")
    '''
    print( 'orig peaks ' + str( xaxis_data[ peaks[0]] ))
    print( 'heights ' + str( peaks[1]['peak_heights']))
    print( 'left peaks ' + str( l_peaks))
    print( 'tx ' + str( tx))
    print( 'lpi ' + str( lpi))
    print( 'lpx ' + str(lpx))
    print( 'lph ' + str(lph))
    '''

    # testing things
    '''
    print( type( xaxis_data))
    print( type( peaks[0]))
    print( peaks[0])
    print( peaks[1]['peak_heights'])
    '''




def only_flight_level_lines( crl_path, crl_name, flight_data_path, flight_name, cutoff_indices, xaxis, xlim='none'):

    warnings.filterwarnings("ignore")

    # update font sizes for this plot!
    TINY_SIZE = 12
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=TINY_SIZE)    # legend fontsize

    # load and process the data
    xr_in_situ = plot_in_situ.load_in_situ( flight_data_path, flight_name, sample_step_size=10)

    # rename variables from xarray for convenience
    str_time = xr_in_situ.str_time
    float_time = xr_in_situ.float_time
    ws = [ float( line) for line in xr_in_situ["WS.d"].values] # convert strings to floats
    wd = [ float( line) for line in xr_in_situ["WD.d"].values]
    uwz = [ float( line) for line in xr_in_situ["UWZ.d"].values]
    rr = [ float( line) for line in xr_in_situ["ASfmrRainRate.1"].values]
    lat = [ float( line) for line in xr_in_situ["LATref"].values ]
    lon = [ float( line) for line in xr_in_situ["LONref"].values ]

    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

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

    else:
        print( 'Error in number of indices! Update them in tc_metadata.py')


    if xaxis == "time":
        xaxis_data = float_time
        xlabel = 'Time (UTC, Hours)'
    elif xaxis == "lon":
        xaxis_data = lon
        xlabel = 'Longitude (Degrees)'
    elif xaxis == "lat":
        xaxis_data = lat
        xlabel = 'Latitude (Degrees)'
    else:
        print( "Please choose 'lon', 'lat', or 'time' as a valid xaxis!")

    # make sure the x axis data is a numpy array for proper plotting
    xaxis_data = np.array(xaxis_data)

    plt.plot( xaxis_data, ws, c='r', label='Tangential Wind Speed (m/s)')
    plt.plot( xaxis_data, uwz, c='k', label='Vertical Wind Speed (m/s)')
    plt.plot( xaxis_data, rr, c='b', label=' SFMR Rain Rate (mm/hr)')

    # ax[0, 0]
    plt.legend(loc='upper left')

    if xlim != 'none':
        plt.xlim( xlim)

    plt.grid('on')

    warnings.filterwarnings("ignore")



def plot_lat_lon( crl_path, crl_name, flight_data_path, flight_name, i1, i2, xaxis, variable_list, xlim='none', title=False):

    # load and process the data
    xr_in_situ = plot_in_situ.load_in_situ( flight_data_path, flight_name, sample_step_size=10)

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
