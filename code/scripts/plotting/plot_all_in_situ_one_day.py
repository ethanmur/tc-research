import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import warnings

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


# main function: cycle through either one tc or all tcs, making crl and in situ plots for all variables
def plot( tc='all'):

    helper_fns.mpl_defaults()

    if tc == 'all':
        tcname_list = [ 'fred', 'grace', 'henri', 'ida', 'sam'] # fred isn't included in all the data bc it doesn't have all the proper metadata!
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        metadata = tc_metadata.all_in_situ_metadata( tcname)

        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['dates'] ))+ '\n')
        # load data

        # initialize current date variable, used for removing duplicate plots
        current_date = metadata['dates'][0]
        current_date_count = 0

        # plot each day's dataset for the tc
        for dataset in range( len( metadata[ 'dates'])):


            # code to skip duplicate cases
            # duplicate case: 2nd statement is to not skip the initial case
            if metadata["dates"][ dataset] == current_date and current_date_count > 0:
                print( 'skipping case ' + str( dataset))
                continue
            # non duplicate case
            else:
                # update metadata to skip next plots
                current_date = metadata['dates'][dataset]
                current_date_count += 1

            # get the correct name of this day's dataset
            crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
            crl_name = tc_metadata.choose_crl_date( metadata['dates'][dataset], metadata['crl_list'])

            title = "All CRL and In Situ Data (Uncorrected) for TC " + metadata['tc_name'] + ", " + metadata[ 'dates'][dataset]


            # make the plots

            # load data: helps with axis trimming
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)

            warnings.filterwarnings("ignore")

            fig = plt.figure( figsize=(27, 13), facecolor='w') # 34, 16

            # plot temperature with original distance axis
            plt.subplot(811)
            plt.title( title)
            make_plots.plot_power_ch1(crl_path, crl_name, index1=0, index2= len( crl_data.time) - 1, xaxis_name='time')

            # plot power with original distance axis
            plt.subplot(812)
            make_plots.plot_T(crl_path, crl_name, index1=0, index2= len( crl_data.time) - 1, xaxis_name='time')

            # get xlims to use later!
            ax = plt.gca()
            ax1lims = ax.get_xlim()


            # add flight level data using newly created in situ netcdf files!
            # load data
            flight_level_path = "/Users/etmu9498/research/data/in-situ-nc"
            flight_level_list = make_plots.load_flight_level( flight_level_path, print_files = False) # all in situ data names
            flight_level_name = tc_metadata.choose_in_situ_date( metadata['dates'][dataset], flight_level_list)
            os.chdir( flight_level_path)
            flight_data = xr.open_dataset( flight_level_name)

            # get the proper axis
            xaxis_data = flight_data.time
            fig = plt.gcf()

            ax1 = fig.add_subplot( 813)
            ax1.set_ylabel('Total Wind Speed (m/s)', c='c')
            ax1.plot( xaxis_data, flight_data['WS.d'], c='c', linewidth=.5)
            ax1.yaxis.tick_right()
            ax1.yaxis.set_label_position("right")
            ax1.grid(True)
            plt.locator_params(axis='x', nbins=10)
            ax1.set_xlim( ax1lims)

            # add an empty colorbar to make everything fit in line
            helper_fns.add_blank_colorbar()

            ax2 = fig.add_subplot( 814)
            degree_symbol = u'\N{DEGREE SIGN}'
            ax2.set_ylabel( 'Wind Dir (' + degree_symbol + ' from ...)')
            ax2.plot( xaxis_data, flight_data['WD.d'], c='g', linewidth=.5)
            ax2.grid(True)
            ax2.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            helper_fns.add_blank_colorbar()

            ax3 = fig.add_subplot( 815)
            ax3.set_ylabel( 'Vertical Wind Speed (m/s)', c='y')
            ax3.plot( xaxis_data, flight_data['UWZ.d'], c='y', linewidth=.5)
            ax3.yaxis.tick_right()
            ax3.yaxis.set_label_position("right")
            ax3.grid(True)
            ax3.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            helper_fns.add_blank_colorbar()

            ax4 = fig.add_subplot( 816)
            ax4.set_ylabel( 'SFMR Rain Rate (mm/hr)', c='b')
            ax4.plot( xaxis_data, flight_data['ASfmrRainRate.1'], c='b', linewidth=.5)
            ax4.grid(True)
            ax4.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            helper_fns.add_blank_colorbar()


            ax4 = fig.add_subplot( 817)
            ax4.set_ylabel( 'P-3 Height (m)', c='r')
            ax4.plot( xaxis_data, flight_data['HT.d'], c='r', linewidth=.5)
            ax4.yaxis.tick_right()
            ax4.yaxis.set_label_position("right")
            ax4.grid(True)

            # find good limits for the final plot
            heights = flight_data['HT.d']

            # new defintion: use surface pressures to choose p-3 height interval?!
            # intense case
            if np.nanmin( flight_data['PSURF.d']) < 850.0:
                minh = 2500
                maxh = np.nanmax( heights[ np.where( heights < 3300.0)]) + 200.0
                continue
            # moderate case
            elif np.nanmin( flight_data['PSURF.d']) < 980.0:
                minh = 2900
                maxh = np.nanmax( heights[ np.where( heights < 3300.0)]) + 50.0
            # weaker case
            elif np.nanmin( flight_data['PSURF.d']) < 1000.0:
                maxh = np.nanmax( heights[ np.where( heights < 3300.0)]) + 25.0
                minh = 3000
            # weak case
            else:
                maxh = np.nanmax( heights[ np.where( heights < 3300.0)]) + 10.0
                minh = 3150

            # the definition below didn't work because the p-3 starts at 0m at takeoff :(
            # so all tcs have a p-3 height of 1000m
            # minh = np.nanmin( heights[ np.where( heights > 1000.0)])
            ax4.set_ylim( [minh, maxh ])# [2900, 3300]) # use 2500 m as a lower limit for strong tcs to see full height dip!

            # print( 'old min: ' + str( np.nanmin( heights[ np.where( heights > 1000.0)])))
            # print( 'new min: ' + str( np.nanmin( heights[ np.where( heights > 2000.0)])))

            ax4.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            helper_fns.add_blank_colorbar()

            # surface pressure data
            ax4 = fig.add_subplot( 818)
            ax4.set_xlabel( 'Time (hours, UTC)')
            ax4.set_ylabel( 'Surf Pressure (hPa)', c='k')
            ax4.plot( xaxis_data, flight_data['PSURF.d'], c='k', linewidth=.5)
            ax4.yaxis.tick_right()
            ax4.yaxis.set_label_position("left")
            ax4.grid(True)

            # find good lower limit for pressure plot
            p = flight_data['PSURF.d']
            maxh = 1023.0 # add some padding to max pressure of 1013.0 hpa
            minh = np.nanmin( p[ np.where( p > 850.0)])
            ax4.set_ylim( [minh - 10.0, maxh])# [2900, 3300]) # use 2500 m as a lower limit for strong tcs to see full height dip!

            ax4.set_xlim(ax1lims)
            plt.locator_params(axis='x', nbins=10)
            helper_fns.add_blank_colorbar()

            print("Plots added")

            os.chdir( "/Users/etmu9498/research/figures/in-situ-all-data")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset) + ".png", dpi=200, bbox_inches='tight')
            print( "Plot " + str( dataset + 1) + " saved\n" )
