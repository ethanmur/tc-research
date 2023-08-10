# import...
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings
import sys
os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/crl-data-processing")
import find_crl_distance_rmws
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops

def plot_all( tc='all', set_xlims=True, save=False, show_limits=False):
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"

    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']

        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( make_plots.load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( make_plots.load_flight_level( crl_data_root + '2022', print_files=False) )

    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2022', print_files=False)]

    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")


    # print out the number of files to be plotted
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))


    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            # use the function below to process the crl data separately!
            crl_data = plot_one( yearval, fileval, set_xlims, save, show_limits)

            print( "New CRL File Plotted and Saved: " +  yearval + "/" + fileval)
    return


def plot_one( yearval, crl_name, set_xlims=True, save=False, show_limits=False):
    # open up the crl file
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"
    crl_path = crl_data_root + yearval
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    plt.figure( figsize=(14, 19))
    helper_fns.change_font_sizes( 12, 12)

    time = crl_data.time
    height = crl_data.height
    
    # plot rmw and radial distance axes!
    plt.subplot( 613)
    plt.plot( time, crl_data.center_dist, c='k', linewidth=1.5)
    plt.ylabel( "Radial Distance (km)")
    ax = plt.gca()
    xlims = ax.get_xlim()
    ax2 = ax.twinx()
    ax2.plot( time, crl_data.rmw, c='g', linewidth=1.5, label='RMW')
    ax2.set_ylabel( "RMW (unitless)")
    # add blank line for colorbar
    ax2.axvline( x = -1000, c='k', label = 'Radial Distance')
    ax2.legend(loc='upper right')
    helper_fns.add_blank_colorbar()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    ax.set_ylim([0, 350])
    ax2.set_ylim([0, 15])
    if xlims[1] > 1:
        plt.xlim( xlims)
    print( "rmw plotted")

    # plot p3 height and temps!
    plt.subplot( 612)
    plt.plot( time, crl_data.p3_height, c='y', linewidth=1.5)
    plt.ylabel( "P-3 Height (m)")
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot( time, crl_data.fl_T, c='r', linewidth=1.5, label='Temp')
    ax2.set_ylabel( "Temperature ( C)")
    ax2.axvline( x = -1000, c='y', label = 'P-3 Height')
    ax2.legend(loc='upper right')
    helper_fns.add_blank_colorbar()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)
    print( "height plotted")

    # plot wind spd and w axes!
    plt.subplot( 611)
    plt.title( "All CRL data for file " + crl_name[:-3])
    plt.plot( time, crl_data.wind_speed, c='c', linewidth=1.5)
    plt.ylabel( "Wind Speed (m/s)")
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot( time, crl_data.w, c='g', linewidth=1.5, label='w')
    ax2.set_ylabel( "W (m/s)")
    # add blank line for colorbar
    ax2.axvline( x = -1000, c='k', label = 'wind speed')
    ax2.legend(loc='upper right')
    helper_fns.add_blank_colorbar()
    # set these x limits no matter what: prevents showing colorbars at -10000!
    if xlims[1] > 1:
        plt.xlim( xlims)
    print( 'speed plotted')

    # plot temperature
    plt.subplot( 615)
    min = 5
    max = 35
    map = plt.cm.get_cmap( "RdYlBu").reversed()
    plt.pcolormesh( time, height, crl_data.T.transpose(), vmin = min, vmax = max, cmap = map)
    plt.colorbar(label="T (Degrees C)")
    plt.ylabel("Height (m)")

    if set_xlims:
        if xlims[1] > 1:
            plt.xlim( xlims)
    ax = plt.gca()
    ax.set_facecolor('k')
    plt.ylim( [ 0, np.nanmax( height)])
    print("temp plotted")

    # plot wv
    plt.subplot( 616)
    min = 0
    max = 25
    plt.pcolormesh( time, height, crl_data.WVMR.transpose(), vmin = min, vmax = max)
    plt.colorbar(label="WVMR (g/kg)")
    plt.ylabel("Height (m)")
    plt.xlabel( "Time (Hours, UTC)")
    if set_xlims:
        if xlims[1] > 1:
            plt.xlim( xlims)
    ax = plt.gca()
    ax.set_facecolor('k')
    plt.ylim( [ 0, np.nanmax( height)])
    print( "wv plotted")


    # plot power ch1
    plt.subplot( 614)
    min = -30
    max = -10
    plt.pcolormesh( time, height, crl_data.P_ch1.transpose(), vmin = min, vmax = max)
    plt.colorbar(label="Power Ch. 1 (dBz)")
    plt.ylabel("Height (m)")

    if set_xlims:
        if xlims[1] > 1:
            plt.xlim( xlims)
    ax = plt.gca()
    ax.set_facecolor('k')
    plt.ylim( [ 0, np.nanmax( height)])
    print( "power plotted")

    # optional: plot eyewall limits here!
    if show_limits:
        date = crl_name[7:11]
        metadata = eyewall_metadata.all_metadata()

        # testing
        # print(yearval)
        # print(crl_name)
        # print(date)
        
        # check if this date exists... if not, give it some empty eyewall limits!
        # also check for fred am and pm cases!!
        if date == '0812':
            if crl_name[11:13] == "H1":
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
            elif crl_name[11:13] == "H2":
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
        elif date in metadata[ yearval]['eyewall_limits'].keys():
            eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
        else:
            eyewall_limits = [ ()]

        print(eyewall_limits)

        eyewall0 = True
        # cycle through eyewall limits and add a red line!
        for eyewall_pair in eyewall_limits:
            for eyelim in eyewall_pair:
                if eyewall0:
                    plt.axvline(x=eyelim, c='r', linewidth=1.25, label="Eye Limit")
                    eyewall0=False
                else:
                    plt.axvline(x=eyelim, c='r', linewidth=1.25)
        plt.legend(loc='upper right')

    if save:
        savedir = "/Users/etmu9498/research/figures/CRL-all-data-processed/" + yearval
        os.chdir( savedir)
        plt.savefig( crl_name[:-3] + ".png", dpi=100, bbox_inches='tight')



# another reasonably simple plotting script originally saved under "code/tests/"
# this script zooms in on individual eye passes and plots them with cloud heights and limits
# spacing is the time in hours to leave data on either side of the eyewall!
def plot_individual_passes(tc='all', lw=2, spacing=.5, save=False):
    # load all 2021 and 2022 metadata here for testing
    metadata = eyewall_metadata.all_metadata()

    # load data and make example plots for all eyewalls!
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"

    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( make_plots.load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( make_plots.load_flight_level( crl_data_root + '2022', print_files=False) )
    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2022', print_files=False)]
    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")

    print("List of Input Files:")
    print(filelist)
        
    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))
    print( 'figures saved to: /Users/etmu9498/research/figures/CRL-all-data-processed/')


    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):    
        for filei, fileval in enumerate( filelist[ yeari]):
            ######
            ## new code: grab the limits for this case
            ######
            date = fileval[7:11]
            
            # check if this date exists... if not, give it some empty eyewall limits!
            # also check for fred am and pm cases!!
            if date == '0812':
                if fileval[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                elif fileval[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
            else:
                eyewall_limits = [ ()]
            
            # make plots for proposed eyewall limits!
            for pairi, pairval in enumerate( eyewall_limits):
           
                print("Year " + yearval + ", " + date + ", pass " + str(pairi))
                # skip empty cases
                if pairval == ():
                    continue
                    
                print( "Eyewall limits: " + str( pairval))
                crl_path = crl_data_root + yearval
                os.chdir( crl_path)
                crl_data = xr.open_dataset( fileval)
                time = crl_data.time
                height = crl_data.height
                        
                # find left and right eyewall locations. Trim the data down to this pair +- a bit extra for nice, faster plotting!
                # only do this if there are eyewall limits to look at
                if len( pairval) > 1:
                    lefteyewall = pairval[ 0]
                    righteyewall = pairval[ 1]
                    li = np.argmin( np.abs( time.values - lefteyewall + spacing))
                    ri = np.argmin( np.abs( time.values - righteyewall - spacing))

                    # make sure the indices don't go over either limit!
                    if li < 0:
                        li = 0
                    if ri > len( time) - 1:
                        ri = len( time) - 1            
                        
                # make figure
                plt.figure( figsize=(11, 9))
                helper_fns.change_font_sizes( 9, 9)

                print("starting figure")
                
                #################
                # find 50 km limits for the eyes of weak TCs! and print them out below
                inds49 = np.where( np.rint( crl_data.center_dist.values) == 49)[0]
                inds50 = np.where( np.rint( crl_data.center_dist.values) == 50)[0]
                inds51 = np.where( np.rint( crl_data.center_dist.values) == 51)[0]
                inds = inds49.tolist() + inds50.tolist() + inds51.tolist()            
                
                
                # plot rmw and radial distance axes!
                plt.subplot( 511)
                plt.title( "All CRL data for file " + fileval[:-3])
                plt.plot( time, crl_data.center_dist, c='k', linewidth=1.5, label='Radial Distance')
                plt.ylabel( "Radial Distance (km)")
                plt.ylim( [-10, 150])
                ax = plt.gca()
                # ax2 = ax.twinx()
                # ax2.plot( time, crl_data.rmw, c='g', linewidth=1.5, label='RMW')
                # ax2.set_ylabel( "RMW (unitless)")
                # plt.ylim( [-.2, 10])
                # add blank line for colorbar
                # ax2.axvline( x = -1000, c='k', label = 'Radial Distance')
                helper_fns.add_blank_colorbar()
                plt.xlim( [ time[li], time[ri]])
                
                # add 50 km time indices as dots!
                timelist = [] # append times here for better, clearer sorting
                for indi, indval in enumerate( inds):
                    # plt.scatter( time[ inds], crl_data.center_dist[ inds], c='g', s=30)
                    plt.axvline( x=time[indval], c='b', linewidth = 1.0)
                    timelist.append( np.round( time[ indval].values, 4))
                # print the list for the user!!
                # print( np.sort( np.array( timelist)).tolist())        
                # add legend info
                plt.axvline( x=-1000, c='b', linewidth = 1.0, label='50 km cutoff')
                plt.axvline( x = lefteyewall, c='g', linewidth=lw, label='chosen limits')
                plt.axvline( x = righteyewall, c='g', linewidth=lw)
                plt.legend(loc='upper right')
                           
                # plot p3 height and temps!
                plt.subplot( 512)
                plt.plot( time, crl_data.p3_height, c='y', linewidth=1.5)
                plt.ylabel( "P-3 Height (m)")
                ax = plt.gca()
                ax2 = ax.twinx()
                ax2.plot( time, crl_data.fl_T, c='r', linewidth=1.5, label='Temp')
                ax2.set_ylabel( "Temperature ( C)")
                ax2.axvline( x = -1000, c='y', label = 'P-3 Height')
                helper_fns.add_blank_colorbar()
                plt.xlim( [ time[li], time[ri]])
                plt.axvline( x = lefteyewall, c='c', linewidth=lw, label='chosen limits')
                plt.axvline( x = righteyewall, c='c', linewidth=lw)
                ax2.legend(loc='upper right')
                # plt.xlim( [ lefteyewall - .5, righteyewall + .5])

                # plot wind spd and w axes!
                plt.subplot( 513)
                plt.plot( time, crl_data.wind_speed, c='c', linewidth=1.5)
                plt.ylabel( "Wind Speed (m/s)")
                ax = plt.gca()
                ax2 = ax.twinx()
                ax2.plot( time, crl_data.w, c='g', linewidth=1.5, label='w')
                ax2.set_ylabel( "W (m/s)")
                # add blank line for colorbar
                ax2.axvline( x = -1000, c='k', label = 'wind speed')
                ax2.legend(loc='upper right')
                helper_fns.add_blank_colorbar()
                plt.xlim( [ time[li], time[ri]])
                print('wind speed plot created')

                # plot temperature
                plt.subplot( 514)
                min = 5
                max = 35
                map = plt.cm.get_cmap( "RdYlBu").reversed()
                plt.pcolormesh( time[li : ri], height, crl_data.T.transpose()[ :, li:ri], vmin = min, vmax = max, cmap = map)
                plt.colorbar(label="T (Degrees C)")
                plt.ylabel("Height (m)")
                ax = plt.gca()
                ax.set_facecolor('k')
                plt.xlim( [ time[li], time[ri]])
                plt.ylim( [ np.nanmin( height), np.nanmax( height)])
                print( "Temperature plot created")
                    
                '''
                # plot wv
                plt.subplot( 615)
                min = 0
                max = 25
                plt.pcolormesh( time[ li:ri], height, crl_data.WVMR.transpose()[ :, li:ri], vmin = min, vmax = max)
                plt.colorbar(label="WVMR (g/kg)")
                plt.ylabel("Height (m)")
                ax = plt.gca()
                ax.set_facecolor('k')
                plt.ylim( [ np.nanmin( height), np.nanmax( height)])
                plt.xlim( [ time[li], time[ri]])                
                print("wv plot created")
                '''
                    
                # plot power ch1
                # also find and plot cloud top heights!!
                if yearval == '2021':
                    min = -30
                elif yearval == '2022':
                    min = -40
                print( "cutoff power = " + str( min))

                H = crl_data.height
                power = crl_data.P_ch1[ li:ri, :]
                axis = crl_data.time[ li:ri]
                p3_height = crl_data.p3_height[ li:ri]
                cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = min)
                plt.subplot( 515)
                max = -10
                plt.pcolormesh( axis, H, power.transpose(), vmin = min, vmax = max)
                plt.plot( cloudtime, cloudheights, c='r', linewidth=1.0, label="cloud tops")
                plt.colorbar(label="Power Ch. 1 (dBz)")
                plt.ylabel("Height (m)")
                ax = plt.gca()
                ax.set_facecolor('k')
                plt.xlabel( "Time (Hours, UTC)")
                # plot proposed eyewall limits!
                plt.ylim( [ np.nanmin( height), np.nanmax( height)])
                plt.xlim( [ time[li], time[ri]])     
                plt.axvline( x = lefteyewall, c='c', linewidth=lw, label='chosen limits')
                plt.axvline( x = righteyewall, c='c', linewidth=lw)
                plt.legend(loc='upper right')   
                print( "power plot created")
                
                if save:
                    savedir = "/Users/etmu9498/research/figures/CRL-all-data-processed/" + yearval + "-eyewall-limits"
                    os.chdir( savedir)
                    plt.savefig( fileval[:-3] + '-' + str( pairi) + ".png", dpi=200, bbox_inches='tight')
                    print( "Figure Saved")                
            print( "File " + str( filei) + " created.\n")
