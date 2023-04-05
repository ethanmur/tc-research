# import...
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings
import sys
import pandas as pd
import datetime
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns


def make_one_plot( year, fl_name, plot_type='all', ri_cond = 30):

    # load the full ships dataset
    print( "Loading SHIPS dataset")
    os.chdir(  "/Users/etmu9498/research/data/ships/")
    file1 = open('lsdiaga_1982_2021_sat_ts_5day.txt', 'r')
    Lines = file1.readlines()
    # there are a bunch of lines in the total ships file!!
    print( "SHIPS dataset created")
    print( "Number of lines in dataset: " + str( len( Lines)))

    fl_path = "/Users/etmu9498/research/data/in-situ-noaa-processed/" + year
    os.chdir( fl_path)
    data = xr.open_dataset( fl_name)

    strname = fl_name[ 11:15].upper()
    # cut off the string at the underscore! allow for 2 and 3 character tc names
    strtemp = ''
    for chari, charval in enumerate( strname):
        if charval == '_':
            break
        else:
            strtemp += charval
    strname = strtemp

    # convenient variables for plotting
    dates = []
    time_since_start = []
    datetimes = []
    date_only = []

    vmax = []
    psurf = []
    sheardir = []
    shearmag = []
    startdate = 0
    starttime = 0

    header_inds = []
    # go through all the lines
    for ind in range( len( Lines)):
        # get the heading lines, and look for this TC's cases!
        if 'HEAD' and strname  in Lines[ ind]:
            # only keep 2021 cases
            if Lines[ ind][ 6 : 8] == str( year[ 2:4]):
                header_inds.append( ind)

    # print valid header indices
    # print( header_inds)
    # do this for all headers
    for headeri, headerval in enumerate( header_inds):
        # add times to the list! increments of 6 hours
        if headeri == 0:
            time_since_start.append( 0)
            # append starting dates and times
            for i in range( headerval,  len( Lines) ):
                if 'HEAD' in Lines[ i]:
                    startdate = Lines[i][6:12]
                    starttime = Lines[i][13:15]
                    break
            print( 'start date and time updated')
        else:
            # otherwise, find the most recent time and add 6 hours!
            time_since_start.append( time_since_start[-1] + 6)

        # add dates
        for i in range( headerval,  len( Lines) ):
            if 'HEAD' in Lines[ i]:
                dates.append( Lines[i][6:12] )
                break
        # add datetime objects!!
        for i in range( headerval,  len( Lines) ):
            if 'HEAD' in Lines[ i]:
                month = int( Lines[i][8:10] )
                day = int( Lines[i][10:12] )
                hours = int( Lines[i][13:15] )
                datetime_orig = datetime.datetime( int( year), month, day, hours)
                datetimes.append( datetime_orig.strftime( "%m/%d %Hh"))
                date_only.append( datetime_orig.strftime( "%m/%d"))
                break
        # search for vmax!
        for i in range( headerval,  len( Lines) ):
            if 'VMAX' in Lines[i]:
                vmax.append( int( Lines[i][12 : 15]) ) # the last 3 vals for vmax at 0 hours
                break
        # repeat for pressure
        for i in range( headerval,  len( Lines) ):
            if 'MSLP' in Lines[i]:
                psurf.append( int( Lines[i][11:15] )  ) # the last 3 vals for vmax at 0 hours
                break
        # repeat for shear mag and direction
        for i in range( headerval,  len( Lines) ):
            if 'SHRD' in Lines[i]:
                shearmag.append( np.round( float( Lines[i][11:15] ) / 10, 2) ) # the last 3 vals for vmax at 0 hours
                break
        for i in range( headerval,  len( Lines) ):
            if 'SHTD' in Lines[i]:
                sheardir.append( int( Lines[i][11 : 15]) ) # the last 3 vals for vmax at 0 hours
                break


    # make the plots!
    if plot_type == 'all':
        plt.figure( figsize = (12, 14))
        helper_fns.change_font_sizes( 14, 14)
        lw = 2
        n = 2  # Keep every 4th label

        plt.subplot( 411)
        plt.title( "SHIPS Derived Plots for TC " + fl_name.title() + ", " + year)
        plt.plot( datetimes, vmax, c='k', linewidth=lw)
        plt.ylabel( "Max Wind Speed (m/s)")
        plt.xticks(rotation=45, ha="right")
        ax = plt.gca()
        fig = plt.gcf()
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
        # old method
        # for label in ax.xaxis.get_ticklabels()[::4]:
        #     label.set_visible(False)
        fig.tight_layout()
        plt.grid()
        #######
        ## 3/7/23 new code
        ## check for RI regions
        #######
        labelcase = True
        ri_inds = [] # save indices for ri onset and ends here!
        for vi in range( 4, len( vmax) ):
            intensification = vmax[ vi] - vmax[ vi - 4]
            if intensification > ri_cond:
                ri_inds.append( (vi-4, vi))

        # go through the original inds and join repeats! since (2, 6) would overlap with (3, 7),
        # change it to (2, 7) for plotting
        # save the joined index pairs below
        ri_inds_join = []
        # loop through these values- if a new join case, overwrite them after saving the index!!
        joinstart = 0
        joinend = 0

        #print( ri_inds)
        for ind in range( len( ri_inds) - 1):

            # base case: set up the indices for the first run
            if ind == 0:
                joinstart = ri_inds[ ind][ 0]
                joinend = ri_inds[ ind][ 1]

            # see if there's any overlap between the end of the current pair and the start of the next index pair
            # if so, keep the starting index the same, but update the last index
            if ri_inds[ ind][ 1] >= ri_inds[ ind+1][0]:
                joinend = ri_inds[ ind+1][1]
                #print( 'join case')

                # last case: add whatever is left to the inds_join list!
                if ind == len( ri_inds) - 2:
                    ri_inds_join.append( ( joinstart, joinend))
                    #print( 'end case')
            # non overlapping case- add the current pair and make the next case the next indices!
            else:
                ri_inds_join.append( ( joinstart, joinend))
                joinstart = ri_inds[ ind+1][0]
                joinend = ri_inds[ ind+1][1]
                #print( 'new case')

            #print( ind)
            #print( (joinstart, joinend))
        #print( ri_inds_join)

        # go through the rapid intesnfication inds and add red ri timeframes to the plots
        for i, indpair in enumerate( ri_inds_join):
            # add a nice label for the legend here! only do this once, though
            if labelcase:
                ax.axvspan( datetimes[ indpair[ 0]], datetimes[ indpair[ 1]], alpha=0.5, color='red', label="RI: +" + str( ri_cond) + "kt in 24h")
                labelcase = False
            else:
                ax.axvspan( datetimes[ indpair[ 0]], datetimes[ indpair[ 1]], alpha=0.5, color='red')
        plt.legend( loc='upper left')


        plt.subplot( 412)
        plt.plot( datetimes, psurf, c='b', linewidth=lw)
        plt.ylabel( "Min Surface Pressure (hPa)")
        plt.xticks(rotation=45, ha="right")
        ax = plt.gca()
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
        plt.grid()

        plt.subplot( 413)
        plt.plot( datetimes, shearmag, c='r', linewidth=lw)
        plt.ylabel( "Shear Mag (m/s)")
        ax = plt.gca()
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
        plt.xticks(rotation=45, ha="right")
        plt.grid()

        plt.subplot( 414)
        plt.plot( datetimes, sheardir, c='y', linewidth=lw)
        plt.ylabel( "Shear Dir (Degrees N)")
        plt.xticks(rotation=45, ha="right")
        plt.grid()
        ax = plt.gca()
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
        # plt.xlabel( "Hours Since " + starttime + " UTC on " + startdate)
        plt.xlabel( "Date ( for " + year + ")")



# plot histograms of ships values! used as a helper function for group_ships_variables and
# can be optionally called
def plot_hists( ships_df, ri_cond = 30):

    data_list = [ ships_df['Vmax (kt)'], ships_df['Intensification (kt)'],
        ships_df['Shear (kt)'], ships_df['Shear Change (kt)'] ]

    subplot_list = [ 221, 222, 223, 224]
    ylabels = [ "Maximum Velocity (kt)", "24h Intensification (kt)", "Shear (kt)", "24h Shear Change (kt)"]
    colors = ['b', 'Grey', 'g', 'y']
    binsize = [ 5, 5, 2, 2]

    plt.figure( figsize=(7.5, 8))
    helper_fns.change_font_sizes(12.5, 12.5)
    plt.suptitle( "Intensity and Shear Distributions from P-3 Flights, 2016-2021")


    for i in range( len( subplot_list)):
        plt.subplot( subplot_list[ i])
        min4hist=np.round(np.nanmin( data_list[ i]),1)-binsize[ i]
        max4hist=np.round(np.nanmax( data_list[ i]),1)+binsize[ i]
        nbins=int((max4hist-min4hist)/binsize[ i])

        plt.hist( data_list[ i],nbins,edgecolor='black', color=colors[ i])
        plt.xlabel( ylabels[ i])
        plt.ylabel('Case Count')

        # add a vertical line representing the ri condition for the intensification subplot
        if subplot_list[ i] == 222:
            plt.axvline( x = ri_cond, c='r', linewidth = 2, label = "RI: +" + str( ri_cond) + "kt in 24h")
            plt.legend()

        fig = plt.gcf()
        fig.tight_layout()
