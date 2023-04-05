# import...
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings
import sys
os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/crl-data-processing")
import find_crl_distance_rmws


def plot_all( tc='all', set_xlims=True):
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"

    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']

        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( make_plots.load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( make_plots.load_flight_level( crl_data_root + '2022', print_files=False) )

    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2022', print_files=False)]

    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")


    # print out the number of files to be plotted
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))


    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            # use the function below to process the crl data separately!
            crl_data = plot_one( yearval, fileval, set_xlims)

            print( "New CRL File Plotted and Saved: " +  yearval + "/" + fileval)
    return


def plot_one( yearval, crl_name, set_xlims=True):
    # open up the crl file
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"
    crl_path = crl_data_root + yearval
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    plt.figure( figsize=(14, 10))
    helper_fns.change_font_sizes( 12, 12)

    time = crl_data.time
    height = crl_data.height

    # plot p3 height and temps!
    plt.subplot( 311)
    plt.title( "All CRL data for file " + crl_name[:-3])
    plt.plot( time, crl_data.p3_height, c='y', linewidth=1.5)
    plt.ylabel( "P-3 Height (m)")
    plt.grid()
    ax = plt.gca()

    xlims = ax.get_xlim()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)

    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot( time, crl_data.fl_T, c='r', linewidth=1.5, label='Temp')
    ax2.set_ylabel( "Temperature ( C)")
    ax2.axvline( x = -1000, c='y', label = 'P-3 Height')
    ax2.legend(loc='upper right')
    plt.grid()
    helper_fns.add_blank_colorbar()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)

    print( "height plotted")

    '''
    # plot rmw and radial distance axes!
    plt.subplot( 411)
    plt.plot( time, crl_data.center_dist, c='k', linewidth=1.5)
    plt.ylabel( "Radial Distance (km)")
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot( time, crl_data.rmw, c='g', linewidth=1.5, label='RMW')
    ax2.set_ylabel( "RMW (unitless)")
    # add blank line for colorbar
    ax2.axvline( x = -1000, c='k', label = 'Radial Distance')
    ax2.legend(loc='upper right')
    plt.grid()
    helper_fns.add_blank_colorbar()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)

    print( "rmw plotted")

    # plot wind spd and w axes!
    plt.subplot( 413)
    plt.plot( time, crl_data.wind_speed, c='c', linewidth=1.5)
    plt.ylabel( "Wind Speed (m/s)")
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot( time, crl_data.w, c='g', linewidth=1.5, label='w')
    ax2.set_ylabel( "W (m/s)")
    # add blank line for colorbar
    ax2.axvline( x = -1000, c='k', label = 'wind speed')
    ax2.legend(loc='upper right')
    plt.grid()
    helper_fns.add_blank_colorbar()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)
    print( 'speed plotted')
    '''

    # plot temperature
    plt.subplot( 312)
    min = 5
    max = 35
    map = plt.cm.get_cmap( "RdYlBu").reversed()
    plt.pcolormesh( time, height, crl_data.T.transpose(), vmin = min, vmax = max, cmap = map)
    plt.colorbar(label="T (Degrees C)")
    plt.ylabel("Height (m)")

    if set_xlims:
        if xlims[1] > 1:
            plt.xlim( xlims)
    ax = plt.gca()
    ax.set_facecolor('k')
    plt.ylim( [ 0, np.nanmax( height)])
    print("temp plotted")

    '''
    # plot wv
    plt.subplot( 615)
    min = 0
    max = 25
    plt.pcolormesh( time, height, crl_data.WVMR.transpose(), vmin = min, vmax = max)
    plt.colorbar(label="WVMR (g/kg)")
    plt.ylabel("Height (m)")
    print( "wv plotted")

    if set_xlims:
        if xlims[1] > 1:
            plt.xlim( xlims)
    ax = plt.gca()
    ax.set_facecolor('k')
    plt.ylim( [ 0, np.nanmax( height)])
    '''

    # plot power ch1

    plt.subplot( 313)
    min = -30
    max = -10
    plt.pcolormesh( time, height, crl_data.P_ch1.transpose(), vmin = min, vmax = max)
    plt.colorbar(label="Power Ch. 1 (dBz)")
    plt.ylabel("Height (m)")

    if set_xlims:
        if xlims[1] > 1:
            plt.xlim( xlims)
    ax = plt.gca()
    ax.set_facecolor('k')
    plt.xlabel( "Time (Hours, UTC)")
    plt.ylim( [ 0, np.nanmax( height)])
    plt.grid('on')
    print( "power plotted")

    savedir = "/Users/etmu9498/research/figures/CRL-all-data-processed/" + yearval
    os.chdir( savedir)
    # plt.savefig( crl_name[:-3] + ".png", dpi=100, bbox_inches='tight')
