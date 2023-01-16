# import functions
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import helper_fns

# this function makes plots highlighting the differences between the two crl distance scripts.
# it looks at all tc cases and saves the results as plots in separate folders
def comparison_plots(tc='all'):
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

        # make a plot for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):

            # load data from new and old sources
            crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])
            inbound_name, outbound_name = tc_metadata.choose_tdr_data( tcname, metadata[ 'tdr_list'], dataset)
            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)
            new_tdr_path = "/Users/etmu9498/research/data/tdr-new"
            new_crl_path = "/Users/etmu9498/research/data/crl-new"
            # load crl data
            os.chdir( metadata[ 'crl_path'])
            crl_data = xr.open_dataset( crl_name)
            os.chdir( new_crl_path)
            new_crl = xr.open_dataset( new_crl_name)
            # make new distance arrays using the new and old methods
            old_dist = find_dist_old_tdr( tcname, dataset)
            new_dist = find_dist_new_tdr( tcname, dataset)

            xtype = metadata[ 'xtype'][dataset]
            if xtype == 'lat':
                xlabels = 'Latitude (Degrees)'
            elif xtype == 'lon':
                xlabels = "Longitude (Degrees)"


            oldi1 = metadata[ 'crl_range'][dataset][0]
            oldi2 = metadata[ 'crl_range'][dataset][1]
            newi1 = 0
            newi2 = len( new_crl.time) - 1

            # comparing lat / lon crl and tdr plots to new distance plots!
            plt.figure( figsize=(18,25))
            color_map = plt.cm.get_cmap( "RdYlBu").reversed()
            helper_fns.change_font_sizes(small=14, medium=14)
            plot_spacing = 30

            # make lat / lon plots for the old dataset first as a baseline
            plt.subplot(611)
            plt.title( 'Lat or Lon Plots Using Inbound and Outbound Datasets')
            make_plots.plot_tdr( metadata['tdr_path'], inbound_name, outbound_name, xtype)
            plt.xlabel( xlabels)
            plt.xlim( metadata[ 'xlims'][ dataset] )
            plt.subplot(612)
            make_plots.plot_T( metadata[ 'crl_path'], crl_name, oldi1, oldi2, xtype)
            plt.xlim( metadata[ 'xlims'][ dataset] )

            # make two plots for crl and tdr distance datasets using old datasets (inbound / outbound)
            plt.subplot( 613)
            plt.title( 'Dist Plots Using Inbound and Outbound Datasets')
            make_plots.plot_tdr( metadata['tdr_path'], inbound_name, outbound_name, 'dist')
            plt.xlabel( "Distance (Km)")
            plt.xlim( [ np.nanmin( old_dist) - plot_spacing, np.nanmax( old_dist) + plot_spacing])
            plt.subplot( 614)
            temp = crl_data.T[oldi1:oldi2, :].where( crl_data.T[oldi1:oldi2, :].values < 50).transpose()
            plt.pcolormesh( old_dist, - crl_data.H, temp, cmap=color_map, vmin=5, vmax=35)
            plt.colorbar()
            ax = plt.gca()
            ax.set_facecolor('k')
            plt.grid('on')
            plt.xlim( [ np.nanmin( old_dist) - plot_spacing, np.nanmax( old_dist) + plot_spacing])

            # make two plots for crl and tdr distance datasets using new tdr dataset
            plt.subplot( 615)
            plt.title( 'Dist Plots Using New Full TDR Dataset')
            make_plots.plot_new_tdr( new_tdr_path, tdr_name, 'dist')
            plt.xlabel( "Distance (Km)")
            plt.xlim( [ np.nanmin( new_dist) - plot_spacing, np.nanmax( new_dist) + plot_spacing])
            plt.subplot( 616)
            temp = new_crl.T[newi1:newi2, :].where( new_crl.T[newi1:newi2, :].values < 50).transpose()
            plt.pcolormesh( new_dist, - new_crl.H, temp, cmap=color_map, vmin=5, vmax=35)
            plt.colorbar()
            ax = plt.gca()
            ax.set_facecolor('k')
            plt.grid('on')
            plt.xlim( [ np.nanmin( new_dist) - plot_spacing, np.nanmax( new_dist) + plot_spacing])

            # save the plots
            os.chdir( "/Users/etmu9498/research/figures/tdr-new-and-old-distance-plots/")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )
            print( "Plot " + str( dataset + 1) + " saved\n" )

    warnings.filterwarnings("default")



# this function uses a longer (but maybe more accurate?) script to define a crl
# distance range. It uses the old inbound / outbound tdr datasets
# it returns an array representing the distance from the center for crl data
def find_dist_old_tdr( tcname, dataset):
    metadata = tc_metadata.all_data( tc= tcname)
    crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])
    inbound_name, outbound_name = tc_metadata.choose_tdr_data( tcname, metadata[ 'tdr_list'], dataset)

    # load crl data
    os.chdir( metadata[ 'crl_path'])
    crl_data = xr.open_dataset( crl_name)
    # load tdr data
    os.chdir( metadata['tdr_path'])
    inbound_data = xr.open_dataset( inbound_name)
    outbound_data = xr.open_dataset( outbound_name)

    # goal 1: find the lat / lon limits for the crl dataset (not eye limits, but just general, non overlap limits!)
    # using lat or lon?
    i1 = metadata[ 'crl_range'][dataset][0]
    i2 = metadata[ 'crl_range'][dataset][1]
    xtype = metadata[ 'xtype'][dataset]
    if xtype == 'lon':
        lim1 = crl_data.Lon[ i1]
        lim2 = crl_data.Lon[ i2]
        inboundx = inbound_data.longitude
        outboundx = outbound_data.longitude
    elif xtype == 'lat':
        lim1 = crl_data.Lat[ i1]
        lim2 = crl_data.Lat[ i2]
        inboundx = inbound_data.latitude
        outboundx = outbound_data.latitude
    else:
        print( 'update the xtype list in tc_metadata.all_data()!')

    # goal 2: try to find the closest lat / lon values in the tdr dataset
    # find the closest values for both inbound and outbound data
    # then, pick the closest dataset from the two!
    i1, x1 = helper_fns.xr_closest_val( inboundx, lim1 )
    i2, x2 = helper_fns.xr_closest_val( outboundx, lim1 )
    # find the closest value between the two!
    # x1 is closest case
    if np.abs( np.subtract( x1, lim1)) < np.abs( np.subtract( x2, lim1)):
        in_i = i1
        in_x = x1
        in_dataset = 'inbound'
    # x2 is closest case
    else:
        in_i = i2
        in_x = x2
        in_dataset = 'outbound'
    # repeat this process for the second crl limit
    i1, x1 = helper_fns.xr_closest_val( inboundx, lim2)
    i2, x2 = helper_fns.xr_closest_val( outboundx, lim2)
    # find the closest value between the two! x1 is closest case
    if np.abs( np.subtract( x1, lim2)) < np.abs( np.subtract( x2, lim2)):
        out_i = i1
        out_x = x1
        out_dataset = 'inbound'
    # x2 is closest case
    else:
        out_i = i2
        out_x = x2
        out_dataset = 'outbound'

    # goal 3: find the distances that correspond to the closest lat / lons found above!
    # the - sign in front of the inbound data is a carry over from make_plots code: one of the axes needed to get flipped!
    # setting the distance for the first crl limit
    if in_dataset == 'inbound':
        in_dist = - inbound_data.radius[ in_i]
    else:
        in_dist = outbound_data.radius[ in_i]
    # repeat for the second crl limit
    if out_dataset == 'inbound':
        out_dist = - inbound_data.radius[ out_i]
    else:
        out_dist = outbound_data.radius[ out_i]
    # goal 4: make a distance array for the crl data given the distance tdr limits
    i1 = metadata[ 'crl_range'][dataset][0]
    i2 = metadata[ 'crl_range'][dataset][1]

    # print( i1)
    # print( i2)
    crl_distance = np.linspace( in_dist, out_dist, num= i2-i1)
    return crl_distance

# this function is shorter and simpler, but I'm having issues with it. It uses the
# new, combined tdr dataset
def find_dist_new_tdr( tcname, dataset):
    # load data
    metadata = tc_metadata.all_data( tc= tcname)
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    tdr_path = "/Users/etmu9498/research/data/tdr-new"
    crl_path = "/Users/etmu9498/research/data/crl-new"
    os.chdir( tdr_path)
    new_tdr = xr.open_dataset( tdr_name)
    os.chdir( crl_path)

    new_crl = xr.open_dataset( crl_name)

    # goal 1: find the lat / lon limits for the crl dataset
    # using lat or lon?
    i1 = 0
    i2 = len( new_crl.time) - 1
    xtype = metadata[ 'xtype'][dataset]
    if xtype == 'lon':
        lim1 = new_crl.Lon[ i1]
        lim2 = new_crl.Lon[ i2]
        tdrx = new_tdr.longitude
    elif xtype == 'lat':
        lim1 = new_crl.Lat[ i1]
        lim2 = new_crl.Lat[ i2]
        tdrx = new_tdr.latitude
    else:
        print( 'update the xtype list in tc_metadata.all_data()!')

    # goal 2: try to find the closest lat / lon values in the tdr dataset
    # this is sooo much easier when only searching through one dataset!
    i1, x1 = helper_fns.xr_closest_val( tdrx, lim1 )
    i2, x2 = helper_fns.xr_closest_val( tdrx, lim2)

    # goal 3: find the distances that correspond to the closest lat / lons found above!
    # setting the distance for the first crl limit
    dist1 = new_tdr.distance[ i1]
    # second limit
    dist2 = new_tdr.distance[ i2]

    # goal 4: make a distance array for the crl data given the distance tdr limits
    # i1 = metadata[ 'crl_range'][dataset][0]
    # i2 = metadata[ 'crl_range'][dataset][1]

    i1 = 0
    i2 = len( new_crl.time) - 1
    # print( i1)
    # print( i2)
    crl_distance = np.linspace( dist1, dist2, num= i2-i1)
    return crl_distance
