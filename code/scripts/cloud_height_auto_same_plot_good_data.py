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
import crl_flight_level_plots


# the auto plot function defaults to
# other valid storm selections include 'all', 'fred', 'grace', 'henri', and 'ida'

# can choose 'tdr' or 'in-situ' for eyewall_alg
def plot( tc='all', eyewall_alg='in-situ'):

    warnings.filterwarnings("ignore")

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    for tcname in tcname_list:

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
            os.chdir( tcdata['tdr_path'] )
            # print( tcdata['tc_name'])
            # print( tcdata[ 'tdr_list'])
            # print( counter)
            inbound_data, outbound_data = tc_metadata.choose_tdr_data( tcdata['tc_name'], tcdata['tdr_list'], counter)

            os.chdir( tcdata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

            os.chdir( tcdata['in_situ_path'] )
            in_situ_name = tc_metadata.choose_in_situ_date( tcdata['dates'][counter], tcdata['in_situ_list'] )

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

            # if the start of the eywall is being determined by tdr data:
            if eyewall_alg == 'tdr':
                # call a function to find the start of the eyewall from the inbound and outbound datasets
                instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis= axis)
            elif eyewall_alg == 'in-situ':
                instartx, outstartx = eyewall_slope.eyewall_start_in_situ( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_name, tcdata['crl_range'][counter], xaxis = axis)
            else:
                print( "Please choose 'tdr' or 'in-situ' for the type of eyewall_alg")
                return

            print( 'Eyewall algorithm complete')
            # add padding to sides of crl and tdr plots to see the data that's being cut off
            step = 1 # .2
            # how thick should the vertical lines be, and what color?
            axvwidth = 4
            axvcolor = 'g'


            plt.figure(figsize=(18, 10), facecolor='w')

            # show tdr and crl case
            if eyewall_alg == 'tdr':
                plt.subplot(211)

                plt.title( "TDR and CRL, TC " + tcdata['tc_name'] + ", "
                    + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

                # plot tdr data
                eyewall_slope_auto.plot_one_tdr_section( tcdata, counter, eyewall_cutoff=False, xaxis = axis, good_data_case=True)

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
                # four index case
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
                plt.axvline( x= instartx, linewidth= axvwidth, c=axvcolor)
                plt.axvline( x= outstartx, linewidth= axvwidth, c=axvcolor)

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

                os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-good-data-tdr/")
                plt.savefig( "tdr-crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )

            # in situ and crl case
            elif eyewall_alg == 'in-situ':

                plt.subplot(211)
                plt.title( "In Situ and CRL, TC " + tcdata['tc_name'] + ", "
                    + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

                plt.axvline( x= instartx, linewidth= axvwidth, c=axvcolor, label='Eyewall Location')
                plt.axvline( x= outstartx, linewidth= axvwidth, c=axvcolor)


                # plot in situ data
                crl_flight_level_plots.only_flight_level_lines( tcdata["crl_path"], crl_data, tcdata['in_situ_path'], in_situ_name, tcdata['crl_range'][counter], axis )
                # plot eyewall cutoffs

                plt.grid('on')

                # mess with limits
                if not np.isnan( instartx) and not np.isnan( outstartx):
                    if instartx < outstartx:
                        plt.xlim( instartx - step, outstartx + step)
                    else:
                        plt.xlim( instartx + step, outstartx - step)

                print( "In Situ Image " + str( counter + 1) + ": complete" )

                plt.subplot( 212)
                # print( "CRL Image " + str( counter + 1) + ": start" )

                # typical case
                if len( tcdata['crl_range'][counter] ) == 2:
                    # plot crl data
                    make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis, show_colorbar=False)
                # four index case
                else:
                    make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis, show_colorbar=False)
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
                    plt.axvline( x= instartx, linewidth= axvwidth, c=axvcolor)
                    plt.axvline( x= outstartx, linewidth= axvwidth, c=axvcolor)


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

                os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-good-data-in-situ/")
                plt.savefig( "in-situ-crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )

            #
            else:
                print( "Error while creating figure in cloud_height_auto_same_plot_good_data.py")

            print( "Plot " + str( counter + 1) + " saved\n" )

            # this code makes a blank figure in between the real figures: it basically
            # just makes the output look prettier
            f,ax = plt.subplots()
            f.set_visible(False)
            f.set_figheight(1) # figure height in inches

    warnings.filterwarnings("default")


# old code to just plot crl data... not super useful so it's been retired here
'''
crl_comparison=False

# plot just crl data case
if not crl_comparison:
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


    # set crl data limits based on eyewall_start() function above
    if not np.isnan( instartx) and not np.isnan( outstartx):
        if instartx < outstartx:
            plt.xlim( instartx - step, outstartx + step)
        else:
            plt.xlim( instartx + step, outstartx - step)

    # plot eyewall cutoffs
    plt.axvline( x= instartx, linewidth= axvwidth, c= axvcolor)
    plt.axvline( x= outstartx, linewidth= axvwidth, c= axvcolor)

    if eyewall_alg == 'tdr':
        os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-good-data/")
        plt.savefig( "hidden-tdr-crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
    elif eyewall_alg == 'in-situ':
        os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-good-data-in-situ/")
        plt.savefig( "hidden-in-situ-crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
    else:
        print( "Error in saving final figure in cloud_height_auto_same_plot_good_data.py")
    print( "Image " + str( counter + 1) + " complete" )

# plotting tdr and crl data case. Note: saves figs to a different location!
else:
'''
