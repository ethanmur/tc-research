import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
import warnings

os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import tc_metadata
import helper_fns
import cloud_height

import sys
sys.path.append("/Users/etmu9498/research/code/scripts/plotting/")
import simple_flight_level_plot





# do some of the same steps as those in the plot() function below, but average the
# w distributions by tc intensity!
# pretty much the same code as above, but with separate distribution profile panels based on intensity
def intensity_figs( width = .1, smoothwidth = 10):
    warnings.filterwarnings("ignore")

    # empty lists that will hold all the height datasets for each intensity category
    td_w, td_cases = [], 0
    ts_w, ts_cases = [], 0
    wh_w, wh_cases = [], 0
    sh_w, sh_cases = [], 0

    # cycle through every dataset
    tcname_list = ['fred', 'grace', 'henri', 'ida', 'sam']
    for tcname in tcname_list:
        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return
        # print some helpful notices to the user
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')

        # get cloud top height information for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):
            # print out the current dataset for the user
            print( "TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset])

            # get current metadata
            crl_path = metadata['um_crl_path']
            fl_path = metadata['new_flight_data_path']
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
            fl_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)
            # load crl data
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)
            # load flight level (fl)
            os.chdir( fl_path)
            fl_data = xr.open_dataset( fl_name)
            # print( 'CRL and fl data loaded')


            # clip crl and fl data to fit eyewall lims
            # find eyewall lims
            distlims = metadata['eyewall_dists_no_eyewalls'][dataset]
            # crl indices
            crl_lim0 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[0] ).values )
            crl_lim1 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[1] ).values )
            # fl indices
            # find the closest fl time to crl time at eye limits
            fl_lim0 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim0] ).values )
            fl_lim1 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim1] ).values )

            # clip down important fields from datasets!
            w = in_situ_to_float( fl_data['UWZ.d'][ fl_lim0 : fl_lim1] )
            # print( 'w distribution found')

            # save w values from this run in the correct intensity category
            cat = metadata['intensity_cat'][dataset] # tc intensity category
            # figure out where to put the height data depending on the tc intensity!
            if cat == 'td':
                td_w += np.ndarray.tolist( w)
                # print( 'td')
                td_cases += 1
            elif cat == 'ts':
                # print( 'ts')
                ts_w += np.ndarray.tolist( w)
                ts_cases += 1
            elif cat == 'wh':
                # print( 'wh')
                wh_w += np.ndarray.tolist( w)
                wh_cases += 1
            elif cat == 'sh':
                # print( 'sh')
                sh_w += np.ndarray.tolist( w)
                sh_cases += 1

    ################
    ## Part 2: making plots
    ################

    # create w distributions for every TC intensity
    # put things in lists for easier looping
    fig_title = [ 'tropical-depressions', 'tropical-storms', 'weak-hurricanes', 'strong-hurricanes']
    # nicely formatted titles for plots, as opposed to saving
    fig_title_nice = [ 'Tropical Depressions', 'Tropical Storms', 'Weak Hurricanes', 'Strong Hurricanes']
    cases = [ td_cases, ts_cases, wh_cases, sh_cases]
    wlist = [ td_w, ts_w, wh_w, sh_w]

    # loop through each case
    for i in range( 4):

        w = wlist[ i]

        # initial variables
        wmin, wmax, wmean = np.nanmin( w), np.nanmax( w), np.nanmean( w)
        w_bin=np.arange( wmin - .5, wmax+.5, width)
        meanw_count = []

        # binning and smoothing w
        # do this for every w bin
        for newi in w_bin:
            # find the points that fall within this height bin for this step
            res=np.where(np.logical_and( w >= newi - width / 2., w <= newi + width / 2. ))
            # append the count to the list!
            meanw_count.append( len( res[0]) / len( w))
        # smooth data before plotting to eliminate noise
        box_pts = smoothwidth
        box = np.ones(box_pts)/box_pts
        prob_smooth = np.convolve( meanw_count, box, mode='same')

        print("Number of cases for " + fig_title_nice[i] + ": " + str(cases[i]))
        plt.clf()
        helper_fns.change_font_sizes(16, 16)
        plt.figure( figsize=(18, 8))

        title = ( "Vertical Velocity Distribution for all " + fig_title_nice[i])
        plt.suptitle( title)

        plt.subplot(121)
        plt.title("Linear Plot")
        plt.ylabel("Probability of a Given W Value (Sums to 100%)")
        plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
        plt.grid('on')

        plt.ylim([-.5, 11])
        plt.xlim([-10, 10])

        # plot the raw and smoothed data
        plt.plot( w_bin, np.array( meanw_count) * 100, linewidth=1, label="Raw Data", c='b', alpha=.35)
        plt.plot( w_bin, prob_smooth * 100, linewidth=2.5, label="Smoothed Data", c='k')
        plt.legend()


        plt.subplot(122)
        plt.title("Log Plot")
        plt.ylabel("Probability of a Given W Value (Sums to 100%)")
        plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
        plt.grid('on')
        plt.yscale('log')

        plt.ylim([.001, 11])
        plt.xlim([-10, 10])


        # plot the raw and smoothed data
        plt.plot(w_bin, np.array( meanw_count) * 100, linewidth=1, label="Raw Data", c='b', alpha=.35)
        plt.plot(w_bin, prob_smooth * 100, linewidth=2.5, label="Smoothed Data", c='k')
        plt.legend()

        # save figure for this intensity!
        os.chdir( "/Users/etmu9498/research/figures/vertical-vels/distributions-intensity")
        plt.savefig( fig_title[ i] + ".png", bbox_inches='tight', dpi=500, transparent=False )


    ########
    ## Part 3: make one last plot with all averaged lines atop one another!
    ########
    # loop through each case again
    plt.clf()
    helper_fns.change_font_sizes(16, 16)
    plt.figure( figsize=(18, 8))

    colors = ['b', 'k', 'y', 'g']
    for i in range( 4):
        w = wlist[ i]
        # initial variables
        wmin, wmax, wmean = np.nanmin( w), np.nanmax( w), np.nanmean( w)
        w_bin=np.arange( wmin - .5, wmax+.5, width)
        meanw_count = []
        # binning and smoothing w
        # do this for every w bin
        for newi in w_bin:
            # find the points that fall within this height bin for this step
            res=np.where(np.logical_and( w >= newi - width / 2., w <= newi + width / 2. ))
            # append the count to the list!
            meanw_count.append( len( res[0]) / len( w))
        # smooth data before plotting to eliminate noise
        box_pts = smoothwidth
        box = np.ones(box_pts)/box_pts
        prob_smooth = np.convolve( meanw_count, box, mode='same')

        # inside loop: add average lines to plot
        plt.subplot(121)
        plt.plot( w_bin, prob_smooth * 100, linewidth=1.5, label=fig_title_nice[i], c=colors[i])
        plt.subplot(122)
        plt.plot( w_bin, prob_smooth * 100, linewidth=1.5, label=fig_title_nice[i], c=colors[i])



    # outside of loop: add plot details and save the figure!
    title = ( "Vertical Velocity Distribution for all TC Intensities")
    plt.suptitle( title)

    plt.subplot(121)
    plt.title("Linear Plot")
    plt.ylabel("Probability of a Given W Value (Sums to 100%)")
    plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
    plt.grid('on')
    plt.legend()

    plt.subplot(122)
    plt.title("Log Plot")
    plt.ylabel("Probability of a Given W Value (Sums to 100%)")
    plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
    plt.grid('on')
    plt.yscale('log')

    # save the final figure!
    os.chdir( "/Users/etmu9498/research/figures/vertical-vels/distributions-intensity")
    plt.savefig( "all-intensities.png", bbox_inches='tight', dpi=500, transparent=False )


















# the same as the function above, but it makes separate intensity plots for inner and outer
# tc w's! It also skips the raw vs smoothed data plots and just shows the distributions.
def intensity_figs_inner_outer( width = .1, smoothwidth = 10):
    warnings.filterwarnings("ignore")

    # empty lists that will hold all the height datasets for each intensity category
    td_w, td_inner, td_outer, td_cases = [], [], [], 0
    ts_w, ts_inner, ts_outer, ts_cases = [], [], [], 0
    wh_w, wh_inner, wh_outer, wh_cases = [], [], [], 0
    sh_w, sh_inner, sh_outer, sh_cases = [], [], [], 0

    # cycle through every dataset
    tcname_list = ['fred', 'grace', 'henri', 'ida', 'sam']
    for tcname in tcname_list:
        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return
        # print some helpful notices to the user
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')

        # get cloud top height information for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):
            # print out the current dataset for the user
            print( "TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset])

            # get current metadata
            crl_path = metadata['um_crl_path']
            fl_path = metadata['new_flight_data_path']
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
            fl_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)
            # load crl data
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)
            # load flight level (fl)
            os.chdir( fl_path)
            fl_data = xr.open_dataset( fl_name)
            # print( 'CRL and fl data loaded')


            # clip crl and fl data to fit eyewall lims
            # find eyewall lims
            distlims = metadata['eyewall_dists_no_eyewalls'][dataset]
            # crl indices
            crl_lim0 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[0] ).values )
            crl_lim1 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[1] ).values )
            # fl indices
            # find the closest fl time to crl time at eye limits
            fl_lim0 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim0] ).values )
            fl_lim1 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim1] ).values )

            # clip down important fields from datasets!
            w = in_situ_to_float( fl_data['UWZ.d'][ fl_lim0 : fl_lim1] )


            ######################
            ## new code: find inner and outer w vels!
            ######################
            # make a radius of max winds, from -1 to 1, to normalize the size of the eye.
            # previous code used the crl limits to do this, but this code uses in situ data...
            # I don't think that's an issue, but maybe check this section again if it is!
            rmwaxis = np.linspace( -1, 1, num=len( w))
            # break the distribution down by rmw distance!
            lowrmw = []
            highrmw = []
            # cycle through every case
            for i in range( len( w)):
                # high rmw case: add cloud heights
                if rmwaxis[ i] > .5 or rmwaxis[ i] < -.5:
                    highrmw.append( w[ i])
                # low rmw case
                else:
                    lowrmw.append( w[ i])

            # save w values from this run in the correct intensity category
            cat = metadata['intensity_cat'][dataset] # tc intensity category
            # figure out where to put the height data depending on the tc intensity!
            if cat == 'td':
                td_w += np.ndarray.tolist( w)
                td_inner += lowrmw
                td_outer += highrmw
                td_cases += 1
            elif cat == 'ts':
                ts_w += np.ndarray.tolist( w)
                ts_inner += lowrmw
                ts_outer += highrmw
                ts_cases += 1
            elif cat == 'wh':
                wh_w += np.ndarray.tolist( w)
                wh_inner += lowrmw
                wh_outer += highrmw
                wh_cases += 1
            elif cat == 'sh':
                sh_w += np.ndarray.tolist( w)
                sh_inner += lowrmw
                sh_outer += highrmw
                sh_cases += 1

    ################
    ## Part 2: making plots
    ################

    # create w distributions for every TC intensity
    # put things in lists for easier looping
    fig_title = [ 'tropical-depressions', 'tropical-storms', 'weak-hurricanes', 'strong-hurricanes']
    # nicely formatted titles for plots, as opposed to saving
    fig_title_nice = [ 'Tropical Depressions', 'Tropical Storms', 'Weak Hurricanes', 'Strong Hurricanes']
    cases = [ td_cases, ts_cases, wh_cases, sh_cases]

    save_names = [ "all-dists.png", 'inner-eye.png', 'outer-eye.png']
    titleend = [ "All Eye Regions", " Inner Eye Regions", ' Outer Eye Regions']

    # do this for total, inner, and outer cases
    for ploti in range( 3):

        # choose the right vertical vels!
        if ploti == 0:
            wlist = [ td_w, ts_w, wh_w, sh_w]
            print( "all eye w data plotted.")
        elif ploti == 1:
            wlist = [ td_inner, ts_inner, wh_inner, sh_inner]
            print( "inner eye w plotted")
        elif ploti == 2:
            wlist = [ td_outer, ts_outer, wh_outer, sh_outer]
            print( "outer eye w plotted")

        # loop through each case again
        plt.clf()
        helper_fns.change_font_sizes(16, 16)
        plt.figure( figsize=(18, 8))

        colors = ['b', 'k', 'y', 'g']
        for i in range( 4):
            w = wlist[ i]
            # initial variables
            wmin, wmax, wmean = np.nanmin( w), np.nanmax( w), np.nanmean( w)
            w_bin=np.arange( wmin - .5, wmax+.5, width)
            meanw_count = []
            # binning and smoothing w
            # do this for every w bin
            for newi in w_bin:
                # find the points that fall within this height bin for this step
                res=np.where(np.logical_and( w >= newi - width / 2., w <= newi + width / 2. ))
                # append the count to the list!
                meanw_count.append( len( res[0]) / len( w))
            # smooth data before plotting to eliminate noise
            box_pts = smoothwidth
            box = np.ones(box_pts)/box_pts
            prob_smooth = np.convolve( meanw_count, box, mode='same')

            # inside loop: add average lines to plot
            plt.subplot(121)
            plt.plot( w_bin, prob_smooth * 100, linewidth=1.5, label=fig_title_nice[i], c=colors[i])
            plt.subplot(122)
            plt.plot( w_bin, prob_smooth * 100, linewidth=1.5, label=fig_title_nice[i], c=colors[i])



        # outside of loop: add plot details and save the figure!
        title = ( "Vertical Velocity Distribution for " + titleend[ ploti])
        plt.suptitle( title)

        plt.subplot(121)
        plt.title("Linear Plot")
        plt.ylabel("Probability of a Given W Value (Sums to 100%)")
        plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
        plt.grid('on')
        plt.legend()

        plt.ylim([-.5, 11])
        plt.xlim([-10, 10])


        plt.subplot(122)
        plt.title("Log Plot")
        plt.ylabel("Probability of a Given W Value (Sums to 100%)")
        plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
        plt.grid('on')
        plt.yscale('log')

        plt.ylim([.001, 11])
        plt.xlim([-10, 10])

        # save the final figure!
        os.chdir( "/Users/etmu9498/research/figures/vertical-vels/distributions-intensity")
        plt.savefig( save_names[ ploti], bbox_inches='tight', dpi=500, transparent=False )








# main function: cycle through either one tc or all tcs, making crl and in situ plots for all variables
def plot( tc='all', average=False, pad=75, width = .1, smoothwidth = 10):

    if tc == 'all':
        tcname_list = [ 'fred', 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        metadata = tc_metadata.all_data( tcname)
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['dates'] ))+ '\n')

        # plot each day's dataset for the tc
        for dataset in range( len( metadata[ 'dates'])):
            # get current metadata
            crl_path = metadata['um_crl_path']
            fl_path = metadata['new_flight_data_path']
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
            fl_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

            # load crl data
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)
            # load flight level (fl)
            os.chdir( fl_path)
            fl_data = xr.open_dataset( fl_name)

            print( 'CRL and fl data loaded')

            # clip crl and fl data to fit eyewall lims
            # find eyewall lims
            distlims = metadata['eyewall_dists_no_eyewalls'][dataset]
            # crl indices
            crl_lim0 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[0] ).values )
            crl_lim1 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[1] ).values )

            # fl indices
            # find the closest fl time to crl time at eye limits
            fl_lim0 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim0] ).values )
            fl_lim1 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim1] ).values )

            # clip down important fields from datasets!
            w = in_situ_to_float( fl_data['UWZ.d'][ fl_lim0 : fl_lim1] )
            fl_dist = fl_data.distance[ fl_lim0 : fl_lim1]
            # used for easier plotting of fl data without having data analysis take forever
            expand_ind = 500
            w_expand = in_situ_to_float( fl_data['UWZ.d'][ fl_lim0 - expand_ind : fl_lim1 + expand_ind] )
            dist_expand = fl_data.distance[ fl_lim0 - expand_ind : fl_lim1 + expand_ind]

            # find eye cloud heights
            cloud_heights, xaxis = cloud_height.find_cloud_heights( crl_name, -30, crl_lim0, crl_lim1, xaxis='in-situ-dist', crl_path=crl_path, new_heights=True)

            print( 'Cloud heights found')

            # make simple comparison plots
            plt.figure( figsize=(15, 14))
            helper_fns.change_font_sizes(16, 16)

            plt.subplot(311)
            title = ( "CRL and Flight Level Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + " Eye Pass " + metadata['eye_pass'][ dataset] )
            plt.title( title)

            make_plots_new_heights.plot_T( crl_path, crl_name, xaxis = 'in-situ-dist')
            plt.xlim( [distlims[0] - pad, distlims[1] + pad])
            plt.axvline( x= distlims[0], c='g', linewidth=4)
            plt.axvline( x= distlims[1], c='g', linewidth=4)

            plt.subplot(312)
            make_plots_new_heights.plot_power_ch1( crl_path, crl_name, xaxis = 'in-situ-dist')

            plt.xlim( [distlims[0] - pad, distlims[1] + pad])
            plt.axvline( x= distlims[0], c='g', linewidth=4)
            plt.axvline( x= distlims[1], c='g', linewidth=4, label='eye limits')

            # also add cloud heights for testing!
            plt.plot( xaxis, cloud_heights, c='r', linewidth = 2.5, label = 'cloud heights')
            plt.legend(loc='upper right')

            plt.subplot(313)
            plt.plot( dist_expand, w_expand, c='y', linewidth=3)

            plt.grid('on')
            plt.xlabel("Distance from TC Center (Km)")
            plt.ylabel("Vertical Velocity (m/s)")
            helper_fns.add_blank_colorbar()

            plt.ylim( [-10, 10])

            plt.xlim( [distlims[0] - pad, distlims[1] + pad])
            plt.axvline( x= distlims[0], c='g', linewidth=4)
            plt.axvline( x= distlims[1], c='g', linewidth=4)

            os.chdir( "/Users/etmu9498/research/figures/vertical-vels/panels")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset) + ".png", dpi=200, bbox_inches='tight')
            print("Panel plot created")


            # manually make a pdf of vertical vel variations!
            plt.clf()
            helper_fns.change_font_sizes(16, 16)

            # initial variables
            wmin, wmax, wmean = np.nanmin( w), np.nanmax( w), np.nanmean( w)
            w_bin=np.arange( wmin - .5, wmax+.5, width)
            meanw_count = []

            # binning and smoothing w
            # do this for every w bin
            for newi in w_bin:
                # find the points that fall within this height bin for this step
                res=np.where(np.logical_and( w >= newi - width / 2., w <= newi + width / 2. ))
                # append the count to the list!
                meanw_count.append( len( res[0]) / len( w))
            # smooth data before plotting to eliminate noise
            box_pts = smoothwidth
            box = np.ones(box_pts)/box_pts
            prob_smooth = np.convolve( meanw_count, box, mode='same')

            # make the plots!
            plt.figure( figsize=(18, 8))

            title = ( "Vertical Velocity Distribution, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + " Eye Pass " + metadata['eye_pass'][ dataset] )
            plt.suptitle( title)

            plt.subplot(121)
            plt.title("Linear Plot")
            plt.ylabel("Probability of a Given W Value (Sums to 100%)")
            plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
            plt.grid('on')

            # plot the raw and smoothed data
            plt.plot(w_bin, np.array( meanw_count) * 100, linewidth=1, label="Raw Data", c='b', alpha=.35)
            plt.plot(w_bin, prob_smooth * 100, linewidth=2.5, label="Smoothed Data", c='k')
            # plt.axvline(x=wmean, c='g', label='mean value')
            plt.legend()


            plt.subplot(122)
            plt.title("Log Plot")
            plt.ylabel("Probability of a Given W Value (Sums to 100%)")
            plt.xlabel("Flight Level Vertical Velocity Values (m/s)")
            plt.grid('on')
            plt.yscale('log')

            # plot the raw and smoothed data
            plt.plot(w_bin, np.array( meanw_count) * 100, linewidth=1, label="Raw Data", c='b', alpha=.35)
            plt.plot(w_bin, prob_smooth * 100, linewidth=2.5, label="Smoothed Data", c='k')
            plt.legend()

            # save the figure
            os.chdir( "/Users/etmu9498/research/figures/vertical-vels/distributions")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset) + ".png", dpi=200, bbox_inches='tight')
            print("Distribution plot created")


            # Make a plot highlighting height and w perturbations
            plt.clf()
            helper_fns.change_font_sizes(16, 16)

            # find and print covariance value for this eye! Also plot height and w perturbs together
            crl_dist = crl_data.in_situ_distance[ crl_lim0:crl_lim1]
            wnew = np.interp( crl_dist, fl_dist, w)
            heights = cloud_heights

            # averaging case: smooth w and h values by the same window and save plots to a different location!
            if average:
                window = 5
                wavg = pd.Series( wnew).rolling(window=window, min_periods=1, center=True).mean()
                havg = pd.Series( cloud_heights).rolling(window = window, min_periods=1, center=True).mean()
                w_height_covar = covar( wavg, havg)
                wprime = simple_perturb( wavg)
                hprime = simple_perturb( havg)
                print( "Cloud height and w covariance = " + str( w_height_covar))

                plt.figure( figsize=(12, 6))
                plt.xlabel( "Distance from TC Center (Km)")
                plt.ylabel( "Perturbation Values")
                plt.grid('on')
                lw=2.5
                plt.plot( crl_dist, wprime, c='k', linewidth=lw, label='w perturbation (m/s)')
                plt.plot( crl_dist, hprime, c='g', linewidth=lw, label='h perturbation (km)')
                plt.legend()
                os.chdir( "/Users/etmu9498/research/figures/vertical-vels/covariances-avg")

            # regular case: just use unsmoothed values
            else:
                w_height_covar = covar( wnew, heights)
                wprime = simple_perturb( wnew)
                hprime = simple_perturb( heights)
                print( "Cloud height and w covariance = " + str( w_height_covar))

                plt.figure( figsize=(12, 6))
                plt.xlabel( "Distance from TC Center (Km)")
                plt.ylabel( "Perturbation Values")
                plt.grid('on')
                lw=2.5
                plt.plot( crl_dist, wprime, c='k', linewidth=lw, label='w perturbation (m/s)')
                plt.plot( crl_dist, hprime, c='g', linewidth=lw, label='h perturbation (km)')
                plt.legend()
                os.chdir( "/Users/etmu9498/research/figures/vertical-vels/covariances")

            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset) + ".png", dpi=200, bbox_inches='tight')
            print( "Plots for dataset " + str( dataset + 1) + " saved\n" )







# exactly the same as the function above, but plots and analysis use time on the
# x axis, not distance, just to double check that nothing is going haywire.
def plot_time( tc='all', average=False, pad=.05):

    if tc == 'all':
        tcname_list = [ 'fred', 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        metadata = tc_metadata.all_data( tcname)
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['dates'] ))+ '\n')

        # plot each day's dataset for the tc
        for dataset in range( len( metadata[ 'dates'])):
            # get current metadata
            crl_path = metadata['um_crl_path']
            fl_path = metadata['new_flight_data_path']
            tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
            fl_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

            # load crl data
            os.chdir( crl_path)
            crl_data = xr.open_dataset( crl_name)
            # load flight level (fl)
            os.chdir( fl_path)
            fl_data = xr.open_dataset( fl_name)

            print( 'CRL and fl data loaded')

            # clip crl and fl data to fit eyewall lims
            # find eyewall lims
            distlims = metadata['eyewall_dists_no_eyewalls'][dataset]
            # crl indices
            crl_lim0 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[0] ).values )
            crl_lim1 = np.argmin( np.abs(crl_data.in_situ_distance - distlims[1] ).values )

            # fl indices
            # find the closest fl time to crl time at eye limits
            fl_lim0 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim0] ).values )
            fl_lim1 = np.argmin( np.abs( fl_data.time - crl_data.time[crl_lim1] ).values )

            # clip down important fields from datasets!
            w = in_situ_to_float( fl_data['UWZ.d'][ fl_lim0 : fl_lim1] )

            fl_dist = fl_data.distance[ fl_lim0 : fl_lim1]
            fl_time = fl_data.time[ fl_lim0 : fl_lim1]

            # used for easier plotting of fl data without having data analysis take forever
            expand_ind = 500
            w_expand = in_situ_to_float( fl_data['UWZ.d'][ fl_lim0 - expand_ind : fl_lim1 + expand_ind] )
            time_expand = fl_data.time[ fl_lim0 - expand_ind : fl_lim1 + expand_ind]

            # find eye cloud heights
            cloud_heights, xaxis = cloud_height.find_cloud_heights( crl_name, -30, crl_lim0, crl_lim1, xaxis='time', crl_path=crl_path, new_heights=True)

            print( 'Cloud heights found')

            # everything above is the same as the first function! time indices were used for
            # comparison, so that big change doesn't even have to be made.


            # make simple comparison plots
            plt.figure( figsize=(15, 14))
            helper_fns.change_font_sizes(16, 16)

            plt.subplot(311)
            title = ( "CRL and Flight Level Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + " Eye Pass " + metadata['eye_pass'][ dataset] )
            plt.title( title)

            make_plots_new_heights.plot_T( crl_path, crl_name, xaxis = 'time')
            plt.xlim( [ crl_data.time[ crl_lim0] - pad, crl_data.time[ crl_lim1] + pad])
            plt.axvline( x= distlims[0], c='g', linewidth=4)
            plt.axvline( x= distlims[1], c='g', linewidth=4)

            plt.subplot(312)
            make_plots_new_heights.plot_power_ch1( crl_path, crl_name, xaxis = 'time')

            plt.xlim( [ crl_data.time[ crl_lim0] - pad, crl_data.time[ crl_lim1] + pad])
            plt.axvline( x= distlims[0], c='g', linewidth=4)
            plt.axvline( x= distlims[1], c='g', linewidth=4, label='eye limits')

            # also add cloud heights for testing!
            plt.plot( xaxis, cloud_heights, c='r', linewidth = 2.5, label = 'cloud heights')
            plt.legend(loc='upper right')

            plt.subplot(313)
            plt.plot( time_expand, w_expand, c='y', linewidth=3)

            plt.grid('on')
            plt.xlabel("Time ( UTC, Hours)")
            plt.ylabel("Vertical Velocity (m/s)")
            helper_fns.add_blank_colorbar()

            plt.ylim( [-10, 10])

            plt.xlim( crl_data.time[ crl_lim0] - pad, crl_data.time[ crl_lim1] + pad)
            plt.axvline( x= crl_data.time[ crl_lim0], c='g', linewidth=4)
            plt.axvline( x= crl_data.time[ crl_lim1], c='g', linewidth=4)

            os.chdir( "/Users/etmu9498/research/figures/vertical-vels/panels-time")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset) + ".png", dpi=200, bbox_inches='tight')
            print("Panel plot created")



            # Make a plot highlighting height and w perturbations
            plt.clf()
            helper_fns.change_font_sizes(16, 16)

            # find and print covariance value for this eye! Also plot height and w perturbs together
            crl_time = crl_data.time[ crl_lim0:crl_lim1]

            # wnew = np.interp( crl_dist, fl_dist, w)
            wnew = np.interp( crl_time, fl_time, w)
            heights = cloud_heights

            # averaging case: smooth w and h values by the same window and save plots to a different location!
            if average:
                window = 5
                wavg = pd.Series( wnew).rolling(window=window, min_periods=1, center=True).mean()
                havg = pd.Series( cloud_heights).rolling(window = window, min_periods=1, center=True).mean()
                w_height_covar = covar( wavg, havg)
                wprime = simple_perturb( wavg)
                hprime = simple_perturb( havg)
                print( "Cloud height and w covariance = " + str( w_height_covar))

                plt.figure( figsize=(12, 6))
                plt.xlabel( "Time (UTC, Hours)")
                plt.ylabel( "Perturbation Values")
                plt.grid('on')
                lw=2.5
                plt.plot( crl_time, wprime, c='k', linewidth=lw, label='w perturbation (m/s)')
                plt.plot( crl_time, hprime, c='g', linewidth=lw, label='h perturbation (km)')
                plt.legend()
                os.chdir( "/Users/etmu9498/research/figures/vertical-vels/covariances-avg-time")

            # regular case: just use unsmoothed values
            else:
                w_height_covar = covar( wnew, heights)
                wprime = simple_perturb( wnew)
                hprime = simple_perturb( heights)
                print( "Cloud height and w covariance = " + str( w_height_covar))

                plt.figure( figsize=(12, 6))
                plt.xlabel( "Time (UTC, Hours)")
                plt.ylabel( "Perturbation Values")
                plt.grid('on')
                lw=2.5
                plt.plot( crl_time, wprime, c='k', linewidth=lw, label='w perturbation (m/s)')
                plt.plot( crl_time, hprime, c='g', linewidth=lw, label='h perturbation (km)')
                plt.legend()
                os.chdir( "/Users/etmu9498/research/figures/vertical-vels/covariances-time")

            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset) + ".png", dpi=200, bbox_inches='tight')
            print( "Plots for dataset " + str( dataset + 1) + " saved\n" )






# helper function used to turn string w values into floats
def in_situ_to_float( return_var):
    return_var_temp = np.zeros( len( return_var))
    for line_ind in range( len( return_var)):
        if return_var[ line_ind] == '':
            return_var_temp[line_ind] = np.nan
        else:
            return_var_temp[ line_ind] = float( return_var[ line_ind])
    return return_var_temp


# this function will hopefully find covariances easily! Taken from Stull 2.4.5a
# var1 and var2 should be lists or numpy arrays, and they must be the same lengths
def covar( var1, var2):
    # make sure variables are numpy arrays
    var1, var2 = np.array( var1), np.array( var2)
    # find the lengths of the variables (N) and mean values
    if len( var1) == len( var2):
        N = len( var1)
    else:
        print( "Covariance error. Please input arrays of equal lengths!")
    mean1, mean2 = np.nanmean( var1), np.nanmean( var2)
    # save
    sum = 0
    for i in range( N - 1):
        sum += ( var1[i] - mean1) * ( var2[i] - mean2)
    return sum / N


# this function is similar to the one above, but it returns an array of perturbations, not a covariance!
def simple_perturb( var):
    # make sure variables are numpy arrays
    var = np.array( var)
    mean = np.mean( var)
    return var - mean
