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
import in_situ_multi_panels
os.chdir(  "/Users/etmu9498/research/code/scripts/save-new-datasets")
import save_crl_data



# pretty much the same function as the one found in ../save-new-datasets/
# testing_distance_coords.py, except it uses distance from tc center values calculated
# from in situ data, not tdr data!
def distance_plots( padding=250, tc='all'):
    warnings.filterwarnings("ignore")

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
        print( "Saving plots to 'in-situ-calculated-distance' folder")

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
            new_crl_path =  "/Users/etmu9498/research/data/crl-new"

            # load crl data
            os.chdir( new_crl_path)
            new_crl = xr.open_dataset( new_crl_name)

            oldi1 = metadata[ 'crl_range'][dataset][0]
            oldi2 = metadata[ 'crl_range'][dataset][1]
            newi1 = 0
            newi2 = len( new_crl.time) - 1

            # accounting for weird Henri case, 8/21, eye 1 where I needed to cut out a flight loop
            four_case = False
            if len( metadata['crl_range'][dataset] ) == 4:
                four_case = True
                oldi3 = metadata['crl_range'][dataset][2]
                oldi4 = metadata['crl_range'][dataset][3]

            color_map = plt.cm.get_cmap( "RdYlBu").reversed()
            helper_fns.change_font_sizes(small=14, medium=14)
            title = ( "Data with New In Situ Derived Distance Scale, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # comparing lat / lon crl and tdr plots to new distance plots!
            plt.figure( figsize=(14, 14), facecolor='w')

            # make plots for crl and tdr distance datasets using new tdr dataset
            plt.subplot( 411)
            plt.title( title)
            make_plots.plot_new_tdr( new_tdr_path, tdr_name, 'dist')
            plt.xlim( [ padding, - padding])

            plt.subplot(412)
            # calculate power backscattered to channel 1
            make_plots.plot_new_power_ch1( new_crl_path, new_crl_name, data_source='in-situ')
            plt.xlim( [padding, - padding])

            plt.subplot( 413)
            make_plots.plot_new_T( new_crl_path, new_crl_name, data_source = 'in-situ')
            plt.xlim( [ padding, - padding])

            new_flight_data_path = metadata['new_flight_data_path']
            new_flight_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

            # plt.subplot( 414)
            in_situ_multi_panels.make_one_subplot(metadata['crl_path'], crl_name, new_flight_data_path, new_flight_name, metadata['crl_range'][dataset])
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

            plt.xlabel( "Distance from TC Center (Km)")
            print( 'Distance Plots Added')

            # save the plots
            os.chdir( "/Users/etmu9498/research/figures/in-situ-calculated-distance/")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=200 )
            print( "Plot " + str( dataset + 1) + " saved\n" )

    warnings.filterwarnings("default")




# the same function as above, but it makes plots using the height corrected CRL data
def distance_plots_new_heights( padding=250, tc='all'):
    warnings.filterwarnings("ignore")

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
        print( "Saving plots to 'in-situ-calculated-distance-new-heights' folder")

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

            new_crl_path = "/Users/etmu9498/research/data/crl-new-matrices"

            # load crl data
            os.chdir( new_crl_path)
            new_crl = xr.open_dataset( new_crl_name)

            newi1 = 0
            newi2 = len( new_crl.time) - 1


            color_map = plt.cm.get_cmap( "RdYlBu").reversed()
            helper_fns.change_font_sizes(small=14, medium=14)
            title = ( "Data with New In Situ Derived Distance Scale, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # comparing lat / lon crl and tdr plots to new distance plots!
            plt.figure( figsize=(18, 18), facecolor='w')

            # make plots for crl and tdr distance datasets using new tdr dataset
            plt.subplot( 411)
            plt.title( title)
            make_plots.plot_new_tdr( new_tdr_path, tdr_name, 'dist')
            plt.xlim( [ padding, - padding])

            plt.subplot(412)
            # calculate power backscattered to channel 1
            make_plots_new_heights.plot_new_power_ch1( new_crl_path, new_crl_name, data_source='in-situ')
            plt.xlim( [padding, - padding])

            plt.subplot( 413)
            make_plots.plot_new_T( new_crl_path, new_crl_name, data_source = 'in-situ')
            plt.xlim( [ padding, - padding])

            new_flight_data_path = metadata['new_flight_data_path']
            new_flight_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)

            plt.subplot( 414)
            in_situ_multi_panels.make_one_subplot(metadata['crl_path'], crl_name, new_flight_data_path, new_flight_name, metadata['crl_range'][dataset])
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

            plt.xlabel( "Distance from TC Center (Km)")
            print( 'Distance Plots Added')

            # save the plots
            os.chdir( "/Users/etmu9498/research/figures/in-situ-calculated-distance/")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( dataset + 1) + " saved\n" )

    warnings.filterwarnings("default")
