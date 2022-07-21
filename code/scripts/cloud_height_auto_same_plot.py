import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import eyewall_slope
import cloud_height
import eyewall_slope_auto

import tc_metadata

# the auto plot function defaults to
# other valid storm selections include 'all', 'fred', 'grace', 'henri', and 'ida'
def plot( tc='sam', tdr_crl=False):
    warnings.filterwarnings("ignore")

    tcdata = tc_metadata.choose_data_eye_cloud_heights( tc)
    if tcdata == 'selected TC name is not yet implemented':
        print( tcdata)
        return

    print( 'number of crl files: ' + str( len( tcdata['xlims'] )))
    for counter in range( len( tcdata[ 'dates'] )):

        if tcdata[ 'xlims'] [ counter][0] > 0.0:
            axis = 'lat'
        else:
            axis = 'lon'

        # load data
        os.chdir( tcdata['tdr_path'] )
        # special case to plot the one crl tdr case that works for ida!
        if tc == 'ida':
            inbound_data = tcdata[ 'tdr_list'] [ 20]
            outbound_data = tcdata['tdr_list'] [ 21]
        else:
            inbound_data = tcdata[ 'tdr_list'] [ counter*2]
            outbound_data = tcdata['tdr_list'] [ counter*2 + 1]
        os.chdir( tcdata['crl_path'] )

        crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

        # code to find cloud heights and plot estimation of top heights
        cutoff_power = -30 # seems to be in the middle of decent power thresholds!
        i1 = tcdata['crl_range'][counter][0]
        i2 = tcdata['crl_range'][counter][1]
        H, time = cloud_height.find_cloud_heights( crl_data, cutoff_power, i1, i2, xaxis=axis)

        # accounting for weird Henri case, 8/21, eye 1 where I needed to cut out a flight loop
        if len( tcdata['crl_range'][counter] ) == 4:
            i3 = tcdata['crl_range'][counter][2]
            i4 = tcdata['crl_range'][counter][3]
            H2, time2 = cloud_height.find_cloud_heights( crl_data, cutoff_power, i3, i4, xaxis=axis)

        # call a function to find the start of the eyewall from the inbound and outbound datasets
        instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis= axis)

        # plot just crl data case
        if not tdr_crl:
            # make a figure and title
            plt.figure(figsize=(18, 5), facecolor='w')
            plt.title( "CRL Data, TC " + tcdata['tc_name'] + ", "
                + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

            # plot crl data
            make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)

            # plot cloud top height data
            plt.scatter( time, H, c= 'r', s=8, marker='o') # s
            plt.plot( time, H, c= 'r', linewidth=1.5, label= 'Cloud Top Height')
            leg = plt.legend(loc='lower left')
            for line in leg.get_lines():
                line.set_linewidth(4.0)

            # set crl data limits based on eyewall_start() function above
            if not np.isnan( instartx) and not np.isnan( outstartx):
                plt.xlim( instartx, outstartx)
            # the original way to set limits: the ones used in choose_data() function!
            # plt.xlim( tcdata['xlims'] [counter][0], tcdata['xlims'] [counter][1])

            os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-auto/" + tcdata['tc_name'].casefold() + "/")
            plt.savefig( "crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
            print( "Image " + str( counter + 1) + " complete" )

        # plotting tdr and crl data case. Note: saves figs to a different location!
        else:
            # print( "TDR Image " + str( counter + 1) + ": start " )
            # make a figure and title
            plt.figure(figsize=(18, 10), facecolor='w')
            plt.subplot(211)
            plt.title( "TDR and CRL, TC " + tcdata['tc_name'] + ", "
                + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

            # plot tdr data
            eyewall_slope_auto.plot_one_tdr_section( tcdata, counter, eyewall_cutoff=False, xaxis = axis)

            step = .1

            if not np.isnan( instartx) and not np.isnan( outstartx):
                if instartx < outstartx:
                    plt.xlim( instartx - step, outstartx + step)
                else:
                    plt.xlim( instartx + step, outstartx - step)
            print( "TDR Image " + str( counter + 1) + ": complete" )

            plt.subplot( 212)
            # print( "CRL Image " + str( counter + 1) + ": start" )

            # typical case
            if len( tcdata['crl_range'][counter] ) == 2:
                # plot crl data
                make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)
            else:
                make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)
                make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][2], tcdata['crl_range'][counter][3], axis, show_colorbar=False)
                plt.scatter( time2, H2, c= 'r', s=8, marker='o') # s
                plt.plot( time2, H2, c= 'r', linewidth=1.5)


            # plot cloud top height data
            plt.scatter( time, H, c= 'r', s=8, marker='o') # s
            plt.plot( time, H, c= 'r', linewidth=1.5, label= 'Cloud Top Height')
            leg = plt.legend(loc='lower left')
            for line in leg.get_lines():
                line.set_linewidth(4.0)

            # plot eyewall cutoffs
            plt.axvline( x= instartx, linewidth=2, c='g')
            plt.axvline( x= outstartx, linewidth=2, c='g')

            # set crl data limits based on eyewall_start() function above
            if not np.isnan( instartx) and not np.isnan( outstartx):
                if instartx < outstartx:
                    plt.xlim( instartx - step, outstartx + step)
                else:
                    plt.xlim( instartx + step, outstartx - step)

            if axis == 'lat':
                plt.xlabel( "Latitude (Degrees)")
            elif axis == 'lon':
                plt.xlabel( "Longitude (Degrees)")

            print( "CRL Image " + str( counter + 1) + ": complete" )

            # the original way to set limits: the ones used in choose_data() function!
            # plt.xlim( tcdata['xlims'] [counter][0], tcdata['xlims'] [counter][1])

            os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-auto/" + tcdata['tc_name'].casefold() + "/")
            plt.savefig( "tdr-crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
            print( "Plot " + str( counter + 1) + " saved" )

        # this code makes a blank figure in between the real figures: it basically
        # just makes the output look prettier
        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(2) # figure height in inches

    warnings.filterwarnings("default")
