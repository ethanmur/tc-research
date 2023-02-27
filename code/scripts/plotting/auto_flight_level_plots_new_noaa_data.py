import xarray as xr
import os
import matplotlib.pyplot as plt
import warnings
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import make_plots
import tc_metadata
import helper_fns
os.chdir("/Users/etmu9498/research/code/scripts/in-situ-scripts")
import in_situ_colorbar_lines
import sys
sys.path.append("/Users/etmu9498/research/code/scripts/plotting/")
import simple_flight_level_plot
import fl_vpeaks_algorithms as peak_algs


# main function: cycle through either one tc or all tcs, making crl and in situ plots for all variables
# inputs:
# tc:
# ylims: a true or false statementdetermining whether or not to use preset y limits when making time series plots.
#        It's nice to set =True when comparing datasets / getting a feel of scale, but =False when looking for
#        error values (anomalously high wvmr, wind speeds, etc)
def plot( tc='all', max_v_requirement=40, window=10, method='pmin_tlim', timelim=30, ylims=True):

    fl_path_root = "/Users/etmu9498/research/data/in-situ-noaa-full/"
    data_path = "/Users/etmu9498/research/data/"

    # make a list of year folders saved in 'in-situ-noaa-full'
    os.chdir( data_path)
    folder_list = [name for name in os.listdir('in-situ-noaa-full')
        if os.path.isdir(os.path.join('in-situ-noaa-full', name))]

    fl_list_total = []
    # go through every year- count the number of datasets!
    for folderi in range( len( folder_list)):
        # go to new directory and get the count + names of files there!
        fl_listi = make_plots.load_flight_level( fl_path_root + folder_list[ folderi], print_files=False)
        fl_list_total += fl_listi
    fl_list_count = len( fl_list_total)


    # use the whole list, or make plots for just one year!
    if tc == 'all':
        print( 'Total Number of plots to be created: ' + str(fl_list_count)+ '\n')
        folder_input = folder_list # use all possible year folders

    # one year case- check the first two characters: a proper year input?
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):

        # check if this year has a folder
        folder_present = False
        for folder in folder_list:
            # the folder exists
            if folder == tc:
                folder_input = [ tc]
                folder_present = True
                fl_listi = make_plots.load_flight_level( fl_path_root + tc, print_files=False)
                print( 'Total Number of plots to be created: ' + str( len( fl_listi))+ '\n')

        # no folder present case:
        if folder_present == False:
            print( "No folder present. Please download flight level data and add a new folder for it!")
            return

    # input names of specific runs as a dict! each year has its own entry
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        folder_input = list( tc.keys())

    # made it out of the loop?
    else:
        print( 'implement the else case! cannot handle individual plots yet.')




    # do this for every year applicable: all years, one year, or a special case!
    for yeari in folder_input:
        fl_path = fl_path_root + yeari
        # normal cases: just load all the datasets from the given year!
        if tc=='all' or ( len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20')):
            fl_list = make_plots.load_flight_level( fl_path, print_files=False)

        # dictionary case: load only the datasets provided!
        elif type( tc) == type( {}):
            fl_list = tc[ yeari]

        # make sure a save folder exists for the output figures!
        # from goes_gifs_2023_update.py
        os.chdir( "/Users/etmu9498/research/figures/in-situ-all-data-new-noaa/time-series")
        if not os.path.isdir( yeari):
            os.makedirs( yeari)
            print( 'New folder created: time-series/' + yeari)

        os.chdir( "/Users/etmu9498/research/figures/in-situ-all-data-new-noaa/time-series-pretty")
        if not os.path.isdir( yeari):
            os.makedirs( yeari)
            print( 'New folder created: time-series-pretty/' + yeari + "-pretty")

        savedir = "/Users/etmu9498/research/figures/in-situ-all-data-new-noaa/time-series/" + yeari
        prettydir = "/Users/etmu9498/research/figures/in-situ-all-data-new-noaa/time-series-pretty/" + yeari

        # load data
        for dataseti in range( len( fl_list)):
            title = "All Flight Level Data for Flight " + fl_list[ dataseti]

            os.chdir( fl_path)
            # don't decode times: instead, use the time indices to create a time axis!
            fl_data = xr.open_dataset( fl_list[ dataseti], decode_times=False)

            # creating the axis
            # this string holds the start and end times. cut down the interval dataset to get the hours, mins, secs!
            interval_str = fl_data.attrs['TimeInterval']
            h = float( interval_str[0:2])
            m = float( interval_str[3:5])
            s = float( interval_str[6:8])
            start_time = h + m / 60 + s / 3600

            # create the time array manually
            # if this is too slow, try using pandas? see stack overflow link below...
            # https://stackoverflow.com/questions/55648630/how-to-decode-the-time-variable-while-using-xarray-to-load-a-netcdf-file
            time = np.empty( ( len( fl_data['Time'])))
            for timei in range( len( fl_data['Time'])):
                # add to time array
                time[ timei] = start_time + timei / 3600



            # useful for plotting later- just in case pmins aren't found below!
            pmins = []
            time_lims = []

            # find where wind speed peaks will be calculated in the figure!
            # smooth wind speeds and surface pressures- more realistic peaks
            window = window # seconds
            spd_avg = pd.Series( fl_data['WS.d']).rolling(window=window, min_periods=1, center=True).mean()
            # test 1: wind speed threshold with no eyewall limits
            # find max TC strength from averaged wind fields- all eye passes from this tc
            vmax = np.nanmax( spd_avg)
            # check based on max TC strength

            peaks = False

            # peaks meet the 50 m/s threshold
            if vmax > max_v_requirement:
                peaks = True
                # find max vel peaks using the min peak height (defined above) for this intensity

                # old method
                # vpeaks = find_peaks_simple( spd_avg)
                # vpeaks = find_peaks_union( spd_avg)
                if method == 'pmin':
                    vpeaks, pmins = peak_algs.find_peaks_pressure_mins( fl_data['PSURF.d'], spd_avg, window)
                elif method == 'pmin_tlim':
                    vpeaks, pmins, time_lims = peak_algs.find_peaks_pmin_time_limit( fl_data['PSURF.d'], spd_avg, time, window, timelim=timelim)
            else:
                peaks = False




            # make figures! yay!
            warnings.filterwarnings("ignore")
            fig = plt.figure( figsize=(24, 17), facecolor='w')

            helper_fns.change_font_sizes( 18, 18)

            ax1 = fig.add_subplot( 811)
            ax1.set_ylabel('Total Wind Speed (m/s)', c='c')
            ax1.plot( time, fl_data['WS.d'], c='c', linewidth=.5)
            ax1.yaxis.tick_right()
            ax1.yaxis.set_label_position("right")
            ax1.grid(True)
            plt.locator_params(axis='x', nbins=10)

            lw1 = .75
            lw4 = 2.0

            # add peak info
            if peaks:
                ax1.set_title( title)
                for peakind in vpeaks:
                    ax1.axvline(x=time[ peakind], c='k', linewidth = lw1)
            else:
                # ax1.text( time[0] + .5, 30, "Wind speeds are too low :(", fontsize=18)
                ax1.set_title( title + " - Wind speeds are too low for analysis")

            if ylims:
                ax1.set_ylim( [0, 80])


            plt.xlim([ time[ 0], time[-1]])
            ax1lims = ax1.get_xlim()
            ax1.set_xlim( ax1lims)
            plt.gca().axes.get_xaxis().set_visible(False)


            ax2 = fig.add_subplot( 812)
            degree_symbol = u'\N{DEGREE SIGN}'
            ax2.set_ylabel( 'Wind Dir (' + degree_symbol + ' from ...)')
            ax2.plot( time, fl_data['WD.d'], c='g', linewidth=.5)
            ax2.grid(True)
            ax2.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            plt.gca().axes.get_xaxis().set_visible(False)
            if ylims:
                ax2.set_ylim( [-5, 365]) # add some spacing to the 0-360 degree limits


            ax3 = fig.add_subplot( 813)
            ax3.set_ylabel( 'Vertical Wind Speed (m/s)', c='y')
            ax3.plot( time, fl_data['UWZ.d'], c='y', linewidth=.5)
            ax3.yaxis.tick_right()
            ax3.yaxis.set_label_position("right")
            ax3.grid(True)
            ax3.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            plt.gca().axes.get_xaxis().set_visible(False)
            if ylims:
                ax3.set_ylim( [-10, 10])

            ax4 = fig.add_subplot( 814)
            ax4.set_ylabel( 'Rain Rate (mm/hr)', c='b')
            ax4.plot( time, fl_data['SfmrRainRate.1'], c='b', linewidth=.5)
            ax4.grid(True)
            ax4.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            plt.gca().axes.get_xaxis().set_visible(False)
            if ylims:
                ax4.set_ylim( [ -3, 60])


            ax5 = fig.add_subplot( 815)
            ax5.set_ylabel( 'Theta E (K)', c='y')
            ax5.plot( time, fl_data['THETAE.d'], c='y', linewidth=.5)
            ax5.grid(True)
            ax5.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            ax5.yaxis.set_label_position("right")
            ax5.yaxis.tick_right()
            plt.gca().axes.get_xaxis().set_visible(False)
            if ylims:
                ax5.set_ylim( [310, 380])

            ax6 = fig.add_subplot( 816)
            ax6.set_ylabel( 'Temp (C)', c='k')
            ax6.plot( time, fl_data['TA.d'], c='k', linewidth=.5)
            ax6.grid(True)
            ax6.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            plt.gca().axes.get_xaxis().set_visible(False)

            ax7 = fig.add_subplot( 817)
            ax7.set_ylabel( 'WVMR (g/kg)', c='r')
            ax7.plot( time, fl_data['MR.d'], c='r', linewidth=.5)
            ax7.grid(True)
            ax7.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            ax7.yaxis.set_label_position("right")
            ax7.yaxis.tick_right()
            plt.gca().axes.get_xaxis().set_visible(False)
            if ylims:
                ax7.set_ylim( [0, 30])



            # surface pressure data
            ax4 = fig.add_subplot( 818)
            ax4.set_xlabel( 'Time (hours, UTC)')
            ax4.set_ylabel( 'Surf Pressure (hPa)', c='k')
            ax4.plot( time, fl_data['PSURF.d'], c='k', linewidth=.5)
            ax4.yaxis.set_label_position("left")
            ax4.grid(True)

            # find good lower limit for pressure plot
            p = fl_data['PSURF.d']
            maxh = 1023.0 # add some padding to max pressure of 1013.0 hpa
            minh = np.nanmin( p[ np.where( p > 850.0)])
            ax4.set_ylim( [minh - 10.0, maxh])# [2900, 3300]) # use 2500 m as a lower limit for strong tcs to see full height dip!

            ax4.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)

            # add peak info
            for peakind in pmins:
                ax4.axvline(x=time[ peakind], c='b', linewidth = lw4)
            for lim in time_lims:
                ax4.axvline(x=time[ lim], c='r', linewidth = lw4)

            # do the same for the wind speed plot
            for peakind in pmins:
                ax1.axvline(x=time[ peakind], c='b', linewidth = lw1)
            for lim in time_lims:
                ax1.axvline(x=time[ lim], c='r', linewidth = lw1)

            # add legends after the fact! use lines way off the figure scales
            ax1.axvline(x=-10, c='k', linewidth=4, label='eyewall limits')
            ax1.axvline(x=-10, c='b', linewidth=4, label='pressure center')
            ax1.axvline(x=-10, c='r', linewidth=4, label='inner core 30 m limit')
            ax1.legend(loc='upper left', fontsize=12.5)
            ax4.axvline(x=-10, c='b', linewidth=4, label='pressure center')
            ax4.axvline(x=-10, c='r', linewidth=4, label='inner core 30 m limit')
            ax4.legend(loc='upper left', fontsize=12.5)

            os.chdir( savedir)
            plt.savefig( fl_list[ dataseti][:-3] + ".png", dpi=200, bbox_inches='tight')
            print( "Plot " + str( dataseti + 1) + " saved\n" )





            # repeat the steps above to save a 'pretty' version of the same figure lol
            fig = plt.figure( figsize=(16, 8), facecolor='w')
            helper_fns.change_font_sizes( 18, 18)

            ax1 = fig.add_subplot( 211)

            ax1.set_ylabel('Total Wind Speed (m/s)', c='c')
            ax1.plot( time, fl_data['WS.d'], c='c', linewidth=1.5)
            ax1.yaxis.tick_right()
            ax1.yaxis.set_label_position("right")
            ax1.grid(True)
            plt.locator_params(axis='x', nbins=10)

            plt.xlim([ time[ 0], time[-1]])
            ax1lims = ax1.get_xlim()
            ax1.set_xlim( ax1lims)

            lw1 = 2.0
            lw4 = 2.0
            # add peak info
            if peaks:
                ax1.set_title( title)
                for peakind in vpeaks:
                    ax1.axvline( x=time[ peakind], c='k', linewidth = lw1)
            else:
                ax1.set_title( title + " - Wind speeds are too low for analysis")

            # surface pressure data
            ax4 = fig.add_subplot( 212)
            ax4.set_xlabel( 'Time (hours, UTC)')
            ax4.set_ylabel( 'Surf Pressure (hPa)', c='k')
            ax4.plot( time, fl_data['PSURF.d'], c='k', linewidth=1.)
            ax4.yaxis.set_label_position("left")
            ax4.grid(True)

            # find good lower limit for pressure plot
            p = fl_data['PSURF.d']
            maxh = 1023.0 # add some padding to max pressure of 1013.0 hpa
            minh = np.nanmin( p[ np.where( p > 850.0)])
            ax4.set_ylim( [minh - 10.0, maxh])# [2900, 3300]) # use 2500 m as a lower limit for strong tcs to see full height dip!

            ax4.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)

            # add peak info
            for peakind in pmins:
                ax4.axvline(x=time[ peakind], c='b', linewidth = lw4)
            for lim in time_lims:
                ax4.axvline(x=time[ lim], c='r', linewidth = lw4)


            # add legends after the fact! use lines way off the figure scales
            ax1.axvline(x=-10, c='k', linewidth=4, label='eyewall limits')
            ax1.legend(loc='upper left', fontsize=12.5)
            ax4.axvline(x=-10, c='b', linewidth=4, label='pressure center')
            ax4.axvline(x=-10, c='r', linewidth=4, label='inner core 30 m limit')
            ax4.legend(loc='upper left', fontsize=12.5)

            os.chdir( prettydir)
            plt.savefig( fl_list[ dataseti][:-3] + ".png", dpi=200, bbox_inches='tight')
            print( "Plot " + str( dataseti + 1) + " saved\n" )
