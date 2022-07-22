import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap

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
def plot( tc='all' ):

    warnings.filterwarnings("ignore")

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        # load data
        tcdata = tc_metadata.choose_data_cloud_tops_good_data( tcname)
        if tcdata == 'selected TC name is not yet implemented':
            print( tcdata)
            return
        # print some helpful notices to the user
        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['xlims'] ))+ '\n')

        # make a plot for every eye pass for this tc
        for counter in range( len( tcdata[ 'dates'] )):
            # choose the correct x axis values depending on the xlims given in the metadata
            # all latitudes are above 0 and lons are below 0, since we are looking at the Northern Atlantic / Gulf of Mexico
            if tcdata[ 'xlims'] [ counter][0] > 0.0:
                axis = 'lat'
            else:
                axis = 'lon'

            # load data
            os.chdir( tcdata['tdr_path'] )
            inbound_data, outbound_data = tc_metadata.choose_tdr_data( tcdata['tc_name'], tcdata['tdr_list'], counter)

            os.chdir( tcdata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

            os.chdir( tcdata['in_situ_path'] )
            in_situ_name = tc_metadata.choose_in_situ_date( tcdata['dates'][counter], tcdata['in_situ_list'] )

            # code to find cloud heights and plot estimation of top heights
            cutoff_power = -30 # seems to be a decent power threshold middle ground!
            i1 = tcdata['crl_range'][counter][0]
            i2 = tcdata['crl_range'][counter][1]
            H, time = cloud_height.find_cloud_heights( crl_data, cutoff_power, i1, i2, xaxis=axis)

            # accounting for weird Henri case, 8/21, eye 1 where I needed to cut out a flight loop
            if len( tcdata['crl_range'][counter] ) == 4:
                i3 = tcdata['crl_range'][counter][2]
                i4 = tcdata['crl_range'][counter][3]
                H2, time2 = cloud_height.find_cloud_heights( crl_data, cutoff_power, i3, i4, xaxis=axis)

            # Find the eywall using tdr data:
            tdr_instartx, tdr_outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis= axis)
            # Find the eyewall using in situ data:
            insitu_instartx, insitu_outstartx = eyewall_slope.eyewall_start_in_situ( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_name, tcdata['crl_range'][counter], xaxis = axis)

            print( 'Eyewall algorithm complete')
            # add padding to sides of crl and tdr plots to see the data that's being cut off
            step = .25 # .2

            # correctly determine the x limits to be used for all plots
            # the in situ vs tdr comparisons are done to make sure the step buffer is applied well to each case!
            # if any of these limits are nans, it'll cause a bunch of issues, so check for that
            if not np.isnan( tdr_instartx) and not np.isnan( tdr_outstartx):
                # case where in situ limits are larger than tdr limits
                if abs( tdr_instartx) < abs( insitu_instartx):

                    # case where inner limits are smaller than outer limits
                    if insitu_instartx < insitu_outstartx:
                        lim1 = insitu_instartx - step
                        lim2 = insitu_outstartx + step
                    else:
                        lim1 = insitu_instartx + step
                        lim2 = insitu_outstartx - step

                # case where tdr limits are larger than in situ limits
                else:
                    # case where inner limits are smaller than outer limits
                    if tdr_instartx < tdr_outstartx:
                        lim1 = tdr_instartx - step
                        lim2 = tdr_outstartx + step
                    else:
                        lim1 = tdr_instartx + step
                        lim2 = tdr_outstartx - step


            # make the full tdr, in situ, and crl plot!!
            plt.figure(figsize=(20, 26), facecolor='w')
            plt.subplot(411)

            plt.title( "TDR Data, TC " + tcdata['tc_name'] + ", "
                + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

            # plot tdr data
            eyewall_slope_auto.plot_one_tdr_section( tcdata, counter, eyewall_cutoff=False, xaxis = axis, good_data_case=True)
            # add eyewall lines
            add_vert_lines(insitu_instartx, insitu_outstartx, tdr_instartx, tdr_outstartx)
            # add a legend explaining what the vertical lines represent
            plt.legend(loc='upper left')

            plt.xlim( [lim1, lim2])
            print( "TDR Image " + str( counter + 1) + ": complete" )

            # plot crl power data
            plt.subplot( 412)
            # typical case
            if len( tcdata['crl_range'][counter] ) == 2:
                # plot crl data
                make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis, show_colorbar=True)
            # four index case
            else:
                make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis, show_colorbar=True)
                make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][2], tcdata['crl_range'][counter][3], axis, show_colorbar=False)
                plt.scatter( time2, H2, c= 'r', s=8, marker='o') # s
                plt.plot( time2, H2, c= 'r', linewidth=1.5)

            plt.title( "Backscattered CRL Power")
            # plot cloud top height data
            plt.scatter( time, H, c= 'r', s=8, marker='o') # s
            plt.plot( time, H, c= 'r', linewidth=1.5, label= 'Cloud Top Height')

            # add eyewall lines
            add_vert_lines(insitu_instartx, insitu_outstartx, tdr_instartx, tdr_outstartx)
            plt.xlim( [lim1, lim2])
            print( "CRL Power Image " + str( counter + 1) + ": complete" )

            # plot crl temperature data
            plt.subplot( 413)
            # typical case
            if len( tcdata['crl_range'][counter] ) == 2:
                # plot crl data
                make_plots.plot_T(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis, show_colorbar=True)
            # four index case
            else:
                make_plots.plot_T(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis, show_colorbar=True)
                make_plots.plot_T(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][2], tcdata['crl_range'][counter][3], axis, show_colorbar=False)

            # add eyewall lines
            plt.title( "CRL Derived Temperature")
            add_vert_lines(insitu_instartx, insitu_outstartx, tdr_instartx, tdr_outstartx)
            plt.xlim( [lim1, lim2])
            print( "CRL Temperature Image " + str( counter + 1) + ": complete" )

            # plot in situ data
            plt.subplot( 414)
            plt.title( "In Situ Measurements")
            crl_flight_level_plots.only_flight_level_lines( tcdata["crl_path"], crl_data, tcdata['in_situ_path'], in_situ_name, tcdata['crl_range'][counter], axis )
            # add eyewall lines
            add_vert_lines(insitu_instartx, insitu_outstartx, tdr_instartx, tdr_outstartx)
            plt.grid('on')

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
            plt.xlim( [lim1, lim2])
            print( "In Situ Image " + str( counter + 1) + ": complete" )


            if axis == 'lat':
                plt.xlabel( "Latitude (Degrees)")
            elif axis == 'lon':
                plt.xlabel( "Longitude (Degrees)")

            os.chdir( "/Users/etmu9498/research/figures/all-data/")
            plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( counter + 1) + " saved\n" )

            # this code makes a blank figure in between the real figures: it basically
            # just makes the output look prettier
            f,ax = plt.subplots()
            f.set_visible(False)
            f.set_figheight(1) # figure height in inches

    warnings.filterwarnings("default")


def add_vert_lines(insitu_instartx, insitu_outstartx, tdr_instartx, tdr_outstartx):

    # how thick should the vertical lines denoting the eyewalls be, and what color?
    axvwidth = 4
    axvcolor = 'g'
    tdr_axvcolor = 'c'
    plt.axvline( x= insitu_instartx, linewidth= axvwidth, c=axvcolor, label='In Situ Eyewall Location')
    plt.axvline( x= insitu_outstartx, linewidth= axvwidth, c=axvcolor)
    plt.axvline( x= tdr_instartx, linewidth= axvwidth, c=tdr_axvcolor, label='TDR Eyewall Location')
    plt.axvline( x= tdr_outstartx, linewidth= axvwidth, c=tdr_axvcolor)
