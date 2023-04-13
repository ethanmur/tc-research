# import...
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings
import sys
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/fl-data-processing")
from save_fl_data import save_one_fl_processed
import find_fl_dists_rmws
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/crl-data-processing")
import find_crl_distance_rmws
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops
import find_cloud_top_stats
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/fl-data-compositing")
import fl_vpeaks_algorithms as peak_algs


# make detailed time series plots of fl and crl data for the specified cases!
# helpful for visualizing data + results
def plot( tc='all', zoom=False):
    lw = 2
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"
    yearlist, filelist = helper_fns_winter2023.get_crl_datasets( tc)

    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))


    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):

            crl_path = crl_data_root + yearval
            os.chdir( crl_path)
            crl_data = xr.open_dataset( fileval)
            time = crl_data.time
            height = crl_data.height

            # find left and right cutoff locations if specified by the user!
            if zoom:
                lefteyewall = zoom[ 0]
                righteyewall = zoom[ 1]
                li = np.argmin( np.abs( time.values - lefteyewall))
                ri = np.argmin( np.abs( time.values - righteyewall))
                # make sure the indices don't go over either limit!
                if li < 0:
                    li = 0
                if ri > len( time) - 1:
                    ri = len( time) - 1

                # redifine variables as smaller / fit to zoom window!
                time = crl_data.time[ li : ri]
                dist = crl_data.center_dist[ li : ri]
                rmw = crl_data.rmw[ li : ri]
                fl_temp = crl_data.fl_T[ li : ri]
                fl_height = crl_data.p3_height[ li : ri]
                spd = crl_data.wind_speed[ li : ri]
                w = crl_data.w[ li : ri]
                crl_temp = crl_data.T.transpose()[ :, li:ri]
                crl_wv = crl_data.WVMR.transpose()[ :, li:ri]
                power = crl_data.P_ch1[ li:ri, :]
            # normal case: do no trimming
            else:
                dist = crl_data.center_dist
                rmw = crl_data.rmw
                fl_temp = crl_data.fl_T
                fl_height = crl_data.p3_height
                spd = crl_data.wind_speed
                w = crl_data.w
                crl_temp = crl_data.T.transpose()
                crl_wv = crl_data.WVMR.transpose()
                power = crl_data.P_ch1

            # make figure
            plt.figure( figsize=(15, 11))
            helper_fns.change_font_sizes( 12, 12)

            # plot rmw and radial distance axes!
            plt.subplot( 611)
            plt.title( "All CRL data for file " + fileval[:-3])
            plt.plot( time, dist, c='k', linewidth=1.5)
            plt.ylabel( "Radial Distance (km)")
            plt.ylim( [-10, 300])
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot( time, rmw, c='g', linewidth=1.5, label='RMW')
            ax2.set_ylabel( "RMW (unitless)")
            plt.ylim( [-.2, 10])
            # add blank line for colorbar
            plt.grid()
            helper_fns.add_blank_colorbar()

            # plot p3 height and temps!
            plt.subplot( 612)
            plt.plot( time, fl_height, c='y', linewidth=1.5)
            plt.ylabel( "P-3 Height (m)")
            plt.grid()
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot( time, fl_temp, c='r', linewidth=1.5, label='Temp')
            ax2.set_ylabel( "Temperature ( C)")
            ax2.legend(loc='upper right')
            plt.grid()
            helper_fns.add_blank_colorbar()

            # plot wind spd and w axes!
            plt.subplot( 613)
            plt.plot( time, spd, c='c', linewidth=1.5)
            plt.ylabel( "Wind Speed (m/s)")
            ax = plt.gca()
            ax2 = ax.twinx()
            ax2.plot( time, w, c='g', linewidth=1.5, label='w')
            ax2.set_ylabel( "W (m/s)")
            # add blank line for colorbar
            ax2.legend(loc='upper right')
            plt.grid()
            helper_fns.add_blank_colorbar()
            print('wind speed plot created')

            # plot temperature
            plt.subplot( 614)
            min = 5
            max = 35
            map = plt.cm.get_cmap( "RdYlBu").reversed()
            plt.pcolormesh( time, height, crl_temp, vmin = min, vmax = max, cmap = map)
            plt.colorbar(label="T (Degrees C)")
            plt.ylabel("Height (m)")
            ax = plt.gca()
            ax.set_facecolor('k')
            plt.ylim( [ np.nanmin( height), np.nanmax( height)])
            print( "Temperature plot created")

            # plot wv
            plt.subplot( 615)
            min = 0
            max = 25
            plt.pcolormesh( time, height, crl_wv, vmin = min, vmax = max)
            plt.colorbar(label="WVMR (g/kg)")
            plt.ylabel("Height (m)")
            ax = plt.gca()
            ax.set_facecolor('k')
            plt.ylim( [ np.nanmin( height), np.nanmax( height)])
            print("wv plot created")


            # plot power ch1
            # also find and plot cloud top heights!!
            if yearval == '2021':
                min = -30
            elif yearval == '2022':
                min = -40

            cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( height, power, time, fl_height, cutoff_power = min)
            plt.subplot( 616)
            max = -10
            plt.pcolormesh( time, height, power.transpose(), vmin = min, vmax = max)
            plt.plot( cloudtime, cloudheights, c='r', linewidth=1.0, label="cloud tops")
            plt.legend(loc='upper right')
            plt.colorbar(label="Power Ch. 1 (dBz)")
            plt.ylabel("Height (m)")
            ax = plt.gca()
            ax.set_facecolor('k')
            plt.xlabel( "Time (Hours, UTC)")
            plt.ylim( [ np.nanmin( height), np.nanmax( height)])


            print( "power plot created")

            print( "Figure for " + str( fileval) + " created.\n")



# plot
def plot_rmw_sections( tc='all', rmw_lim='default', save_new_data=False, save_fig=False, trimmed_data=False):
    lw = 2

    if trimmed_data:
        fl_data_root = "/Users/etmu9498/research/data/in-situ-noaa-trimmed/"
    else:
        fl_data_root = "/Users/etmu9498/research/data/in-situ-noaa-processed/"

    yearlist, filelist = helper_fns_winter2023.get_fl_datasets( tc)

    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))


    # do this for all the datasets specified! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):

            fl_path = fl_data_root + yearval
            os.chdir( fl_path)
            fl_data = xr.open_dataset( fileval)

            # split up the dataset into distinct radial passes!
            # for separate plotting

            #########
            # step 1: find velocity peaks for creating the rmw axes for all passes
            # kinda a redundant process based on the code in "2023-02-15 flight tracks to radial dist v3.ipynb"
            # probs a better solution, but this works!
            #########
            max_v_requirement=40
            window=10
            timelim = 60
            spd_avg = pd.Series( fl_data['WS.d']).rolling(window=window, min_periods=1, center=True).mean()
            # find the indices of the velocity peaks, pressure minima, and
            vpeaks, pmins, time_lims = peak_algs.find_peaks_pmin_time_limit( fl_data['PSURF.d'], spd_avg, fl_data.time.values, window, timelim=timelim)

            # edge case: no vpeaks are present / no rmw axis created
            # just save the whole dataset again! no splits
            if len( vpeaks) == 0:
                all_inds = [ 0, len( fl_data.time.values) - 1]
            # normal case: vpeaks present
            else:
                dists = fl_data.center_dist.values
                rmw0_vals = []
                for peak in pmins:
                    rmw0_vals.append( dists[ peak])
                rmw1_vals = []
                for peak in vpeaks:
                    rmw1_vals.append( dists[ peak])

                # divide up the distance axis into distinct RMW regions
                rainbandlist = []
                for vpeaki, vpeak in enumerate( vpeaks):
                    #print( 'vpeak index = ' + str( vpeaki))
                    # figure out the distances to the peak winds and current pressure center
                    rmw1 = rmw1_vals[ vpeaki]
                    # use // to get rid of remainder
                    pressurei = vpeaki // 2
                    rmw0 = rmw0_vals[ pressurei]

                    # first value case: count everything before the first wind speed to the center
                    if vpeaki == 0:
                        rainband_ind = 0
                        rainbandlist.append( rainband_ind)
                        #print( 'rainband index = ' + str( rainband_ind))
                        continue
                    # last value case: count everything from the last center to the end of the dataset
                    elif vpeaki == len( vpeaks) - 1:
                        rainband_ind = len( fl_data.time.values) - 1
                        rainbandlist.append( rainband_ind)
                        #print( 'rainband index = ' + str( rainband_ind))
                        continue
                    # middle cases: find the rainband indices manually!
                    # odd cases: the center has already been passed
                    elif vpeaki % 2 == 1:
                        # find the index in the rainbands
                        # should be halfway between the two closest pressure centers
                        # first pressure center + ( difference / 2)
                        rainband_ind = int( pmins[ pressurei] + ( ( pmins[ pressurei + 1] - pmins[ pressurei]) / 2) )
                        rainbandlist.append( rainband_ind)
                        #print( 'rainband index = ' + str( rainband_ind))
                    # even cases: the center is being approached
                    elif vpeaki % 2 == 0:
                        # find the index in the rainbands
                        rainband_ind = int( pmins[ pressurei - 1] + ( ( pmins[ pressurei ] - pmins[ pressurei - 1]) / 2) )
                        rainbandlist.append( rainband_ind)
                        #print( 'rainband index = ' + str( rainband_ind))

                ########
                # step 2: stick pressure minima indices into the data!
                ########
                ############
                ## 3/17 new code: find the actual smallest rmws without finding vpeaks again!
                # use the scipy.find_peaks() fn
                # height > -.5 ensures that pmins are really low, aka less than .1 RMW.
                # distance > 300 ensures that pmins are more than 300 data points apart...
                # see sam 9/29 for when this is an issue, 2 peaks really close to each other!
                # this ensures that pmins aren't slightly off: the arrays start and end
                # at their smallest points!
                ############
                rmw = fl_data.rmw.values
                pminsnew = find_peaks( - rmw, height=-.5, distance=300)[0]

                #print( pmins)
                #print( pminsnew)
                #print( rainbandlist)

                # combine lists and remove duplicate rainfall inds
                all_inds = list( set( rainbandlist + pminsnew.tolist())) # pmins.tolist()
                # sort in ascending order
                all_inds = list( np.sort( np.array( all_inds)))
                #print( all_inds)

            ##############
            # step 3: go through each ind pair and make a simple test rmw plot!
            ##############

            # add 1d relative vorticity to the dataset / plot!
            # using the method found in kossin and eastin 2001
            # give the helper fn total wind speed in m/s and radial distance in meters
            # return processed / smoothed relative vorticity in units of 10^-4 s-1
            rel_vort = find_fl_dists_rmws.calc_rel_vort( fl_data['WS.d'].values, fl_data['center_dist'].values * 1000, double_average=False, window=60)

            # add it to the dataframe for nicer loops!
            fl_data[ 'rel_vort'] = rel_vort

            # list all variables to plot here
            vars = ['TA.d', 'MR.d', 'WS.d', 'HUM_REL.d', 'rel_vort']
            labels = ["Temperature (C)", "Water Vapor (g/kg)", "Wind Spped (m/s)", "Relative Humidity (%)", "Relative Vorticity (10^-4 s-1)"]
            colors = ['r', 'b', 'c', 'k', 'g']
            ylims = [ [0, 30], [0, 25], [0, 80], [0, 100], [-100, 200]]
            # make figure before loop
            fig = plt.figure( figsize=(6 * len( vars), 3 * len( all_inds)))
            helper_fns.change_font_sizes( 14, 14)
            gs=GridSpec( len( all_inds), len( vars))
            plt.suptitle( "Case " + fileval)
            # fig.supxlabel( "Radius of Maximum Winds")
            # fig.tight_layout(rect=[0, 0.03, 1, 0.95])


            for currenti, currentval in enumerate( all_inds):
                # edge case: stop here!
                if currenti == len( all_inds) - 1:
                    break
                else:
                    # special case: save data!
                    # use the helper function in save_fl_data.py to make the dataset
                    if save_new_data:
                        new_fl = save_one_fl_processed( fl_data, (all_inds[currenti], all_inds[currenti+1]), fileval )
                        # final goal: save the newly created crl dataset
                        # make a name for the dataset
                        filename = fileval[:-3] + "_case" + str( currenti) + ".nc"
                        # make sure the folder exists!
                        os.chdir("/Users/etmu9498/research/data/in-situ-noaa-individual")
                        output_folder = yearval
                        if not os.path.isdir( output_folder):
                            os.makedirs( output_folder)
                            print( 'New folder created: ' + output_folder)
                        # go to the new folder and save the data!
                        new_fl.to_netcdf('/Users/etmu9498/research/data/in-situ-noaa-individual/' + yearval + "/" + filename)
                        print( "New In Situ File Created and Saved: " + filename + "\n")

                    for vari, varval in enumerate( vars):
                        rmw = fl_data.rmw[ all_inds[currenti] : all_inds[currenti+1]]
                        temp = fl_data[ varval][ all_inds[currenti] : all_inds[currenti+1]]
                        # trim out large rmws. calc temp first before changing rmw!
                        temp = temp[ np.where( rmw <= 5)[0]]
                        rmw = rmw[ np.where( rmw <= 5)[0]]

                        plt.subplot(gs[ currenti, vari])
                        plt.plot( rmw, temp, c=colors[ vari], label= "Pass " + str(currenti + 1))
                        plt.legend( loc='upper right')
                        plt.xlim( [0, 5])
                        if ylims[ vari]:
                            plt.ylim( ylims[ vari])
                        if currenti == 0:
                            plt.title( labels[ vari])
                        elif currenti == len( all_inds) - 2:
                            plt.xlabel( "Radius of Maximum Winds")

            if save_fig:
                # make sure a save folder exists for the output figures!
                # from goes_gifs_2023_update.py

                if trimmed_data:
                    path = "/Users/etmu9498/research/figures/in-situ-all-data-new-noaa/time-series-trimmed/"
                else:
                    path = "/Users/etmu9498/research/figures/in-situ-all-data-new-noaa/time-series-individual/"
                
                os.chdir( path)
                if not os.path.isdir( yearval):
                    os.makedirs( yearval)
                    print( 'New folder created: time-series/' + yearval)
                savedir = path + yearval
                os.chdir( savedir)
                plt.savefig( fileval[:-3] + ".png", dpi=200, bbox_inches='tight')
                print( "Plot " + fileval + " saved\n" )
