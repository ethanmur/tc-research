# import functions
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from matplotlib.colors import ListedColormap
import xarray as xr
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import helper_fns
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts")
import in_situ_colorbar_lines
os.chdir(  "/Users/etmu9498/research/code/scripts/save-new-datasets")
import save_crl_data

# this function makes plots highlighting the differences between the two crl distance scripts.
# it looks at all tc cases and saves the results as plots in separate folders
def comparison_plots(tc='all', shift_crl_dist=True):
    warnings.filterwarnings("ignore")

    # rerun the crl saving function to apply new shifts and stretches!
    if shift_crl_dist == True:
        save_crl_data.save_all_crl( shift_crl_dist=True, add_dist_coords=True)

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        # load data
        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return
        # print some helpful notices to the user
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')

        # make a plot for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):

            # load data from new and old sources
            crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])
            inbound_name, outbound_name = tc_metadata.choose_tdr_data( tcname, metadata[ 'tdr_list'], dataset)
            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)
            new_tdr_path = "/Users/etmu9498/research/data/tdr-new"
            new_crl_path = "/Users/etmu9498/research/data/crl-new"
            # load crl data
            os.chdir( new_crl_path)
            new_crl = xr.open_dataset( new_crl_name)
            # make a new distance array using the new method
            new_dist = new_crl.distance # find_dist_new_tdr( tcname, dataset)

            xtype = metadata[ 'xtype'][dataset]
            if xtype == 'lat':
                xlabels = 'Latitude (Degrees)'
            elif xtype == 'lon':
                xlabels = "Longitude (Degrees)"

            oldi1 = metadata[ 'crl_range'][dataset][0]
            oldi2 = metadata[ 'crl_range'][dataset][1]
            newi1 = 0
            newi2 = len( new_crl.time) - 1
            # newi2 = metadata[ 'crl_range'][dataset][1] - metadata[ 'crl_range'][dataset][0] - 1

            # accounting for weird Henri case, 8/21, eye 1 where I needed to cut out a flight loop
            four_case = False
            if len( metadata['crl_range'][dataset] ) == 4:
                four_case = True
                oldi3 = metadata['crl_range'][dataset][2]
                oldi4 = metadata['crl_range'][dataset][3]


            # comparing lat / lon crl and tdr plots to new distance plots!
            plt.figure( figsize=(30,18))
            color_map = plt.cm.get_cmap( "RdYlBu").reversed()
            helper_fns.change_font_sizes(small=14, medium=14)
            # add padding on either end of the x axis (in km)
            plot_spacing = 0
            title = ( "TDR Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )
            plt.suptitle( title)

            # padding for distance plots
            crlx = 100

            # special x lims for lat / lon method:
            lim0 = metadata[ 'xlims'][ dataset][0]
            lim1 = metadata[ 'xlims'][ dataset][1]
            old_lim_scale = 1.0

            # lat case
            if lim0 > 0:
                diff = lim0 - lim1
                midpoint = lim1 + diff / 2

                lim0 = midpoint - old_lim_scale
                lim1 = midpoint + old_lim_scale
            # lon case
            else:
                diff = lim1 - lim0
                midpoint = lim0 + diff / 2

                lim0 = midpoint - old_lim_scale
                lim1 = midpoint + old_lim_scale

            old_lims = [ lim0, lim1]

            # make lat / lon plots for the old dataset first as a baseline
            plt.subplot(321)
            plt.title( 'Lat or Lon Plots Using Inbound and Outbound Datasets')
            make_plots.plot_tdr( metadata['tdr_path'], inbound_name, outbound_name, xtype)
            plt.xlim( old_lims )

            plt.subplot(323)
            make_plots.plot_power_ch1( metadata[ 'crl_path'], crl_name, oldi1, oldi2, xtype)
            if four_case:
                make_plots.plot_power_ch1( metadata[ 'crl_path'], crl_name, oldi3, oldi4, xtype, show_colorbar=False)
            plt.xlim( old_lims )

            plt.subplot(325)
            make_plots.plot_T( metadata[ 'crl_path'], crl_name, oldi1, oldi2, xtype)
            if four_case:
                make_plots.plot_T( metadata[ 'crl_path'], crl_name, oldi3, oldi4, xtype, show_colorbar=False)
            plt.xlim( old_lims)
            plt.xlabel( xlabels)

            print( 'Lat Lon Plots Added')

            # make plots for crl and tdr distance datasets using new tdr dataset
            plt.subplot( 322)
            plt.title( 'Dist Plots Using New Full TDR Dataset')
            make_plots.plot_new_tdr( new_tdr_path, tdr_name, 'dist')
            plt.xlabel( "Distance (Km)")
            # plt.xlim( [ np.nanmin( new_dist) - plot_spacing, np.nanmax( new_dist) + plot_spacing])
            plt.xlim( [ crlx, - crlx])

            plt.subplot(324)
            # calculate power backscattered to channel 1
            make_plots.plot_new_power_ch1( new_crl_path, new_crl_name)
            # plt.xlim( [ np.nanmin( new_dist) - plot_spacing, np.nanmax( new_dist) + plot_spacing])
            plt.xlim( [crlx, - crlx])

            plt.subplot( 326)
            make_plots.plot_new_T( new_crl_path, new_crl_name)
            # plt.xlim( [ np.nanmin( new_dist) - plot_spacing, np.nanmax( new_dist) + plot_spacing])
            plt.xlim( [ crlx, - crlx])
            plt.xlabel( "Distance (Km)")

            print( 'Distance Plots Added')

            # save the plots
            os.chdir( "/Users/etmu9498/research/figures/tdr-new-and-old-distance-plots/")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( dataset + 1) + " saved\n" )

    warnings.filterwarnings("default")




# pretty much the same function as above, but it only saves the distance plots
def distance_plots( padding=250, tc='all', shift_crl_dist=True, in_situ=False):
    warnings.filterwarnings("ignore")

    # rerun the crl saving function to apply new shifts and stretches!
    if shift_crl_dist == True:
        save_crl_data.save_all_crl( shift_crl_dist=True, add_dist_coords=True)

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        # load data
        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return
        # print some helpful notices to the user
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')
        if in_situ:
            print( "Saving plots to 'distance-plots-in-situ' folder")
        elif padding >= 250:
            print( "Saving pltos to 'distance-plots-zoom-out' folder")
        else:
            print( "Saving plots to 'distance-plots-zoom-in' folder")

        # make a plot for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):

            # load data from new and old sources
            crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])
            inbound_name, outbound_name = tc_metadata.choose_tdr_data( tcname, metadata[ 'tdr_list'], dataset)
            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)

            flight_data_path = metadata['in_situ_path']
            flight_data_list = make_plots.load_flight_level( flight_data_path, print_files=False)
            flight_name = tc_metadata.choose_in_situ_date( metadata['dates'][dataset], flight_data_list)

            new_tdr_path = "/Users/etmu9498/research/data/tdr-new"
            new_crl_path = "/Users/etmu9498/research/data/crl-new"
            # load crl data
            os.chdir( new_crl_path)
            new_crl = xr.open_dataset( new_crl_name)
            # make a new distance array using the new method
            new_dist = new_crl.tdr_distance # find_dist_new_tdr( tcname, dataset)

            xtype = metadata[ 'xtype'][dataset]
            if xtype == 'lat':
                xlabels = 'Latitude (Degrees)'
            elif xtype == 'lon':
                xlabels = "Longitude (Degrees)"

            oldi1 = metadata[ 'crl_range'][dataset][0]
            oldi2 = metadata[ 'crl_range'][dataset][1]
            newi1 = 0
            newi2 = len( new_crl.time) - 1
            # newi2 = metadata[ 'crl_range'][dataset][1] - metadata[ 'crl_range'][dataset][0] - 1

            # accounting for weird Henri case, 8/21, eye 1 where I needed to cut out a flight loop
            four_case = False
            if len( metadata['crl_range'][dataset] ) == 4:
                four_case = True
                oldi3 = metadata['crl_range'][dataset][2]
                oldi4 = metadata['crl_range'][dataset][3]

            color_map = plt.cm.get_cmap( "RdYlBu").reversed()
            helper_fns.change_font_sizes(small=14, medium=14)
            # add padding on either end of the x axis (in km)
            title = ( "CRL and TDR Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # adding a fourth in situ plot
            if in_situ:
                # comparing lat / lon crl and tdr plots to new distance plots!
                plt.figure( figsize=(18, 18), facecolor='w')

                # make plots for crl and tdr distance datasets using new tdr dataset
                plt.subplot( 411)
                plt.title( title)
                make_plots.plot_new_tdr( new_tdr_path, tdr_name, 'dist')
                plt.xlim( [ padding, - padding])

                plt.subplot(412)
                # calculate power backscattered to channel 1
                make_plots.plot_new_power_ch1( new_crl_path, new_crl_name, data_source='tdr')
                plt.xlim( [padding, - padding])

                plt.subplot( 413)
                make_plots.plot_new_T( new_crl_path, new_crl_name, data_source = 'tdr')
                plt.xlim( [ padding, - padding])

                plt.subplot( 414)
                cutoff_indices = (newi1, newi2)
                in_situ_colorbar_lines.only_flight_level_lines( new_crl_path, new_crl_name, flight_data_path, flight_name, cutoff_indices, 'dist', tcname=tcname, dataset=dataset)
                plt.xlim( [ padding, - padding])

                # add an empty colorbar to make everything fit in line... kinda a
                # messy solution but it's ok for now!
                viridis = cm.get_cmap('viridis', 256)
                newcolors = viridis(np.linspace(0, 1, 256))
                white = np.array([ 1, 1, 1, 1])
                newcolors[:, :] = white
                white_cmap = ListedColormap(newcolors)

                map = matplotlib.cm.ScalarMappable(cmap= white_cmap, norm=matplotlib.colors.Normalize( vmin= 0, vmax= 1))
                cbar = plt.colorbar( mappable= map)
                cbar.set_ticks([])
                cbar.outline.set_visible(False)
            # only adding three plots
            else:
                # comparing lat / lon crl and tdr plots to new distance plots!
                plt.figure( figsize=(18, 18), facecolor='w')

                # make plots for crl and tdr distance datasets using new tdr dataset
                plt.subplot( 411)
                plt.title( title)
                make_plots.plot_new_tdr( new_tdr_path, tdr_name, 'dist')
                plt.xlim( [ padding, - padding])

                plt.subplot(412)
                # calculate power backscattered to channel 1
                make_plots.plot_new_power_ch1( new_crl_path, new_crl_name)
                plt.xlim( [padding, - padding])

                plt.subplot( 413)
                make_plots.plot_new_T( new_crl_path, new_crl_name)
                plt.xlim( [ padding, - padding])

            plt.xlabel( "Distance from TC Center (Km)")
            print( 'Distance Plots Added')

            # save the plots
            if in_situ:
                os.chdir( "/Users/etmu9498/research/figures/distance-plots-in-situ/")
            # zoom out case
            elif padding >= 250:
                os.chdir( "/Users/etmu9498/research/figures/distance-plots-zoom-out")
            # zoom in case
            else:
                os.chdir( "/Users/etmu9498/research/figures/distance-plots-zoom-in/")

            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( dataset + 1) + " saved\n" )

    warnings.filterwarnings("default")
