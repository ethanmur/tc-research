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


# this function is very similar to plot_rmw.py, except that it doesn't use the rmw
# axis at all! instead, in the second plot, the values chosen for the eyewalls divide
# the distance axis shown in the first plot; the result is that rmw = 1 falls where the
# chosen distances are
def plot( tc='all', padding=50):

    helper_fns.mpl_defaults()


    if tc == 'all':
        tcname_list = [ 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # run this now for nice font sizes later!
    helper_fns.change_font_sizes()

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        metadata = tc_metadata.all_data( tcname)

        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['dates'] ))+ '\n')
        # load data

        # plot each day's dataset for the tc
        for dataset in range( len( metadata[ 'dates'])):
            # get the correct name of this day's dataset
            crl_path = "/Users/etmu9498/research/data/crl-new-matrices"
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)

            title = " RMW Axis Creation, CRL and In Situ Data, TC " + metadata['tc_name'] + ", " + metadata[ 'dates'][dataset] + ' Eye Pass ' + metadata[ 'eye_pass'][dataset]

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

            # get eyewall distances for this case
            dists = metadata['in_situ_eyewall_dists'][ dataset]


            # issue with two positive eyewall distances: can't properly split up the data!
            if dists[0] < 0.0:
                # split distance dataset at 0 km
                # this code is very similar to that in find_crl_distance.py
                dist_pos = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance > 0.0)].values
                dist_neg = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance <= 0.0)].values

                # divide by the limits to find rmw values!
                rmw_pos = dist_pos / dists[1]
                rmw_neg = - dist_neg / dists[0]

                '''
                print( len( rmw_pos))
                print( len( rmw_neg))
                print( len( dist_pos))
                print( len( dist_neg))
                print( len( rmw))
                '''
                # combine datasets again
                rmw = np.concatenate( (rmw_neg, rmw_pos))

                # make xlims correspond: find rmw values closest to padding index
                idx1 = (np.abs(crl_data.in_situ_distance + padding )).argmin().values
                idx2 = (np.abs(crl_data.in_situ_distance - padding )).argmin().values
                rmw_lim1 = rmw[ idx1]
                rmw_lim2 = rmw[idx2]

                # this wasn't working so well for sam, so I switched to a static rmw limit
                # the data doesn't match up as well for some cases now, but it works better
                # for sam
                # plt.xlim( [ rmw_lim1, rmw_lim2])
                plt.xlim( [-2, 2])

                # use these definitions below to properly map the rmw axis onto in situ data!
                rmw_lim_minus1 = ( np.abs( rmw - 1.0)).argmin() # np.where( crl_data.rmw_negatives == -1.0)
                rmw_lim0 = ( np.abs( rmw)).argmin() # np.where( crl_data.rmw_negatives == 0.0)
                rmw_lim_plus1 = ( np.abs( rmw + 1.0)).argmin() # np.where( crl_data.rmw_negatives == 1.0)

                '''
                print( rmw[ idx1])
                print( rmw[ idx2])
                print( crl_data.in_situ_distance[ rmw_lim_minus1].values)
                print( crl_data.in_situ_distance[ rmw_lim0].values)
                print( crl_data.in_situ_distance[ rmw_lim_plus1].values)

                np.set_printoptions(threshold=np.inf)
                # print( rmw)
                np.set_printoptions(threshold=1000)
                '''

                # plot rmw!
                plt.pcolormesh(  rmw, - crl_data.H_new, crl_data.power_new.transpose(), vmin = -30, vmax =-10)
            else:
                print( "Height error: both eyewall distance values are positive :(")

            plt.ylabel( 'Height (km)')
            plt.colorbar(label="Backscattered Power ( dBz)")
            plt.grid( 'on')
            plt.ylim( [ 0, crl_data.H_max + .1])
            plt.locator_params(axis='x', nbins=11)
            plt.xlabel( "Radius of Max Winds (Unitless)")
            ax = plt.gca()
            ax.set_facecolor('k')


            # plot in situ data
            # plt.subplot(313)
            # in_situ_colorbar_lines.only_flight_level_lines( crl_path, crl_name, in_situ_path, in_situ_name, cutoff_indices, 'dist', tcname=tcname, dataset=dataset)
            simple_flight_level_plot.plot( crl_path, crl_name, metadata, dataset, 'dist')
            plt.xlim( [ - padding, padding])
            plt.locator_params(axis='x', nbins=11)
            plt.xlabel( "Distance from TC Center (km)")

            # add lines denoting rmw = 0 and = 1 on in situ data to check peaks:
            # is it just a data spacing error or a actual value errors?
            # statistical calcs would be ok even with spacing / plotting errors, but
            # if the scales are actually off, that would be bad

            '''
            # these lines of code are incorrect! I'm pulling RMW data from the crl dataset, when I really want to be
            # using the rmw axis defined above

            # indices of 3 rmw values of interest
            rmw_lim_minus1 = ( np.abs( crl_data.rmw_negatives - 1.0)).argmin() # np.where( crl_data.rmw_negatives == -1.0)
            rmw_lim0 = ( np.abs( crl_data.rmw_negatives)).argmin() # np.where( crl_data.rmw_negatives == 0.0)
            rmw_lim_plus1 = ( np.abs( crl_data.rmw_negatives + 1.0)).argmin() # np.where( crl_data.rmw_negatives == 1.0)
            '''

            # load in situ data
            in_situ_path = metadata[ 'new_flight_data_path']
            in_situ_name = tc_metadata.choose_new_in_situ_name( metadata['tc_name'], dataset)
            os.chdir( in_situ_path)
            in_situ_data = xr.open_dataset( in_situ_name)


            # get indices of in situ data corresponding to crl data
            idx1 = (np.abs(in_situ_data.distance - crl_data.in_situ_distance[rmw_lim_minus1].values )).argmin().values
            idx2 = (np.abs(in_situ_data.distance - crl_data.in_situ_distance[rmw_lim0].values )).argmin().values
            idx3 = (np.abs(in_situ_data.distance - crl_data.in_situ_distance[rmw_lim_plus1].values )).argmin().values

            ''' looots of testing! to basically track down the bug improperly place in situ axes
            print( in_situ_name)
            print( dists)
            print( 'rmw minus: ' + str( rmw[ rmw_lim_plus1]))
            print( 'rmw zero: ' + str( rmw[ rmw_lim0]))
            print( 'rmw plus: ' + str( rmw[ rmw_lim_minus1]))

            print( 'crl minus: ' + str( crl_data.in_situ_distance[rmw_lim_plus1].values))
            print( 'crl zero: ' + str( crl_data.in_situ_distance[rmw_lim0].values))
            print( 'crl plus: ' + str( crl_data.in_situ_distance[rmw_lim_minus1].values))

            print( 'in situ minus: ' + str( in_situ_data.distance[ idx3].values))
            print( 'in situ zero: ' + str( in_situ_data.distance[ idx2].values))
            print( 'in situ plus: ' + str( in_situ_data.distance[ idx1].values))
            '''
            
            # add vert lines showing peak
            plt.axvline( x = in_situ_data.distance[ idx1], c='g', alpha=.75)
            plt.axvline( x = in_situ_data.distance[ idx2], c='g', alpha=.75)
            plt.axvline( x = in_situ_data.distance[ idx3], c='g', alpha=.75)

            # remove this line
            helper_fns.add_blank_colorbar()

            print( "Plot " + str( dataset + 1) + " created\n" )

            os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/rmw-chosen-eyewalls")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=200 )
            print( "Plot " + str( dataset + 1) + " saved\n" )









# The same as the function above, except that it only compares outputs from
# chosen rmws vs calculated rmws!
def plot_differences( tc='all', padding=50):

    helper_fns.mpl_defaults()


    if tc == 'all':
        tcname_list = [ 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # run this now for nice font sizes later!
    helper_fns.change_font_sizes()

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        metadata = tc_metadata.all_data( tcname)

        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['dates'] ))+ '\n')
        # load data

        # plot each day's dataset for the tc
        for dataset in range( len( metadata[ 'dates'])):
            # get the correct name of this day's dataset
            crl_path = "/Users/etmu9498/research/data/crl-new-matrices"
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)

            title = " RMW Axis Creation, CRL and In Situ Data, TC " + metadata['tc_name'] + ", " + metadata[ 'dates'][dataset] + ' Eye Pass ' + metadata[ 'eye_pass'][dataset]

            # make the plots

            # load data: helps with axis trimming
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)

            # get eyewall distances for this case
            dists = metadata['in_situ_eyewall_dists'][ dataset]

            warnings.filterwarnings("ignore")

            fig = plt.figure( figsize=(20, 16), facecolor='w')

            # plot power with original distance axis
            plt.subplot(311)
            plt.title( title)
            make_plots_new_heights.plot_power_ch1(crl_path, crl_name, 'in-situ-dist')
            plt.xlim( [ - padding, padding])
            plt.xlabel( "Distance from TC Center (km)")
            plt.ylim( [ 0, crl_data.H_max + .1])

            # two vertical lines showing where eyewall lims are on distance plot
            plt.axvline( x= dists[0], c='r', alpha=.75, linewidth=3)
            plt.axvline( x= dists[1], c='r', alpha=.75, linewidth=3)

            # plot power with new hand chosen rmw axis
            plt.subplot(312)

            # issue with two positive eyewall distances: can't properly split up the data!
            if dists[0] < 0.0:
                # split distance dataset at 0 km
                # this code is very similar to that in find_crl_distance.py
                dist_pos = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance > 0.0)].values
                dist_neg = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance <= 0.0)].values

                # divide by the limits to find rmw values!
                rmw_pos = dist_pos / dists[1]
                rmw_neg = - dist_neg / dists[0]

                # combine datasets again
                rmw = np.concatenate( (rmw_neg, rmw_pos))

                # make xlims correspond: find rmw values closest to padding index
                idx1 = (np.abs(crl_data.in_situ_distance + padding )).argmin().values
                idx2 = (np.abs(crl_data.in_situ_distance - padding )).argmin().values
                rmw_lim1 = rmw[ idx1]
                rmw_lim2 = rmw[idx2]
                plt.xlim( [ rmw_lim1, rmw_lim2])
                rmw_lim_minus1 = ( np.abs( rmw - 1.0)).argmin() # np.where( crl_data.rmw_negatives == -1.0)
                rmw_lim0 = ( np.abs( rmw)).argmin() # np.where( crl_data.rmw_negatives == 0.0)
                rmw_lim_plus1 = ( np.abs( rmw + 1.0)).argmin() # np.where( crl_data.rmw_negatives == 1.0)

                # plot rmw!
                plt.pcolormesh(  rmw, - crl_data.H_new, crl_data.power_new.transpose(), vmin = -30, vmax =-10)
            else:
                print( "Height error: both eyewall distance values are positive :(")

            plt.ylabel( 'Height (km)')
            plt.colorbar(label="Backscattered Power ( dBz)")
            plt.grid( 'on')
            plt.ylim( [ 0, crl_data.H_max + .1])
            plt.locator_params(axis='x', nbins=11)
            plt.xlabel( "Radius of Max Winds ( From Hand Picked Eyewalls)")
            ax = plt.gca()
            ax.set_facecolor('k')


            # plot calculated, not chosen, rmw axis
            plt.subplot(313)
            make_plots_new_heights.plot_power_ch1( crl_path, crl_name, 'rmw-negatives')
            plt.ylim( [ 0, crl_data.H_max + .1])
            # make xlims correspond: find rmw values closest to padding index
            idx1 = (np.abs(crl_data.in_situ_distance + padding )).argmin().values
            idx2 = (np.abs(crl_data.in_situ_distance - padding )).argmin().values

            rmw_lim1 = crl_data.rmw_negatives[ idx1]
            rmw_lim2 = crl_data.rmw_negatives[idx2]
            plt.locator_params(axis='x', nbins=11)
            plt.xlim( [ rmw_lim1, rmw_lim2])
            plt.xlabel( "Radius of Max Winds ( From In Situ Max WS)")



            print( "Plot " + str( dataset + 1) + " created\n" )

            os.chdir( "/Users/etmu9498/research/figures/CRL-new-heights/rmw-comparisons")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=200 )
            print( "Plot " + str( dataset + 1) + " saved\n" )
