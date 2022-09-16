import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import warnings

os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import tc_metadata
import helper_fns
os.chdir("/Users/etmu9498/research/code/scripts/in-situ-scripts")
import in_situ_colorbar_lines
import sys
sys.path.append("/Users/etmu9498/research/code/scripts/plotting/")
import simple_flight_level_plot


# this function tests the rmw axes by creating 3 subplots: the first two compare
# in situ data with distance and rmw axes, while the third plots in situ data to
# show where the max winds are. The axes limits are all scaled so that each figure
# is spatially proportional. Vertical lines on the in situ plot denote distances where
# rmw = -1, 0, and 1
def rmw_test( tc='all', padding=50):

    helper_fns.mpl_defaults()


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

            title = " RMW Axis Creation, CRL and In Situ Data, TC " + tcdata['tc_name'] + ", " + tcdata[ 'dates'][counter] + ' Eye Pass ' + tcdata[ 'eye_pass'][counter]

            # make the plots

            # load data: helps with axis trimming
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)

            warnings.filterwarnings("ignore")

            fig = plt.figure( figsize=(20, 14), facecolor='w')

            # plot power with original distance axis
            plt.subplot(311)
            plt.title( title)
            make_plots_new_heights.plot_power_ch1(crl_path, crl_name, 'in-situ-dist')
            plt.xlim( [ - padding, padding])
            plt.xlabel( "Distance from TC Center (km)")
            plt.ylim( [ 0, crl_data.H_max + .1])

            # plot power with new rmw axis
            plt.subplot(312)
            make_plots_new_heights.plot_power_ch1( crl_path, crl_name, 'rmw-negatives')

            plt.ylim( [ 0, crl_data.H_max + .1])
            # make xlims correspond: find rmw values closest to padding index
            idx1 = (np.abs(crl_data.in_situ_distance + padding )).argmin().values
            idx2 = (np.abs(crl_data.in_situ_distance - padding )).argmin().values

            rmw_lim1 = crl_data.rmw_negatives[ idx1]
            rmw_lim2 = crl_data.rmw_negatives[idx2]
            plt.locator_params(axis='x', nbins=11)
            plt.xlim( [ rmw_lim1, rmw_lim2])
            plt.xlabel( "Radius of Max Winds (Unitless)")

            # plot in situ data
            # plt.subplot(313)
            # in_situ_colorbar_lines.only_flight_level_lines( crl_path, crl_name, in_situ_path, in_situ_name, cutoff_indices, 'dist', tcname=tcname, dataset=counter)
            simple_flight_level_plot.plot( crl_path, crl_name, tcdata, counter, 'dist')
            plt.xlim( [ - padding, padding])
            plt.locator_params(axis='x', nbins=11)
            plt.xlabel( "Distance from TC Center (km)")

            # add lines denoting rmw = 0 and = 1 on in situ data to check peaks:
            # is it just a data spacing error or a actual value errors?
            # statistical calcs would be ok even with spacing / plotting errors, but
            # if the scales are actually off, that would be bad

            # indices of 3 rmw values of interest
            rmw_lim_minus1 = ( np.abs( crl_data.rmw_negatives - 1.0)).argmin() # np.where( crl_data.rmw_negatives == -1.0)
            rmw_lim0 = ( np.abs( crl_data.rmw_negatives)).argmin() # np.where( crl_data.rmw_negatives == 0.0)
            rmw_lim_plus1 = ( np.abs( crl_data.rmw_negatives + 1.0)).argmin() # np.where( crl_data.rmw_negatives == 1.0)

            # load in situ data
            in_situ_path = tcdata[ 'new_flight_data_path']
            in_situ_name = tc_metadata.choose_new_in_situ_name( tcdata['tc_name'], counter)
            os.chdir( in_situ_path)
            in_situ_data = xr.open_dataset( in_situ_name)
            # get indices of in situ data corresponding to crl data
            idx1 = (np.abs(in_situ_data.distance - crl_data.in_situ_distance[rmw_lim_minus1].values )).argmin().values
            idx2 = (np.abs(in_situ_data.distance - crl_data.in_situ_distance[rmw_lim0].values )).argmin().values
            idx3 = (np.abs(in_situ_data.distance - crl_data.in_situ_distance[rmw_lim_plus1].values )).argmin().values

            # add vert lines showing peak
            plt.axvline( x = in_situ_data.distance[ idx1], c='g', alpha=.75)
            plt.axvline( x = in_situ_data.distance[ idx2], c='g', alpha=.75)
            plt.axvline( x = in_situ_data.distance[ idx3], c='g', alpha=.75)

            # remove this line
            helper_fns.add_blank_colorbar()

            print( "Plot " + str( counter + 1) + " created\n" )

            os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/rmw-tests")
            plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png", bbox_inches='tight', dpi=150 )
            print( "Plot " + str( counter + 1) + " saved\n" )
