# Author: Ethan Murray
# Created:
# Edited: 3/8/23
# Helper functions for save_fl_data.py

# import...
import numpy as np
import os
import sys
import pandas as pd
from geopy import distance
import metpy.calc as mpcalc
os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/fl-data-compositing")
import fl_vpeaks_algorithms as peak_algs



# calculate 1d relative vorticity for flight level data using the definition found in kossin and eastin 2001!
# vtan is the total wind speed (m/s) and r is the radial distance from the center (m).
# double averaging makes the profiles look a lot more like those found in the paper, but idk why really...
# window = 20 seems to be a good middle ground window
def calc_rel_vort( vtan, r, double_average=False, window=20):
    rel_vort = []
    # find relative vorticity for every velocity point! except the last to avoid i+1 error
    for i in range( len( vtan) - 1):
        vi = vtan[ i] # current vtan
        vi1 = vtan[ i+1] # the next vtan value
        vavg = (vi + vi1) / 2
        ri = r[ i] # current radial distance
        ri1 = r[ i+1] # next radial distance
        ravg = (ri + ri1) / 2
        rel_vort.append( vavg / ravg + (vi1 - vi)/(ri1 - ri) )
    # add one last nan value to get arrays to be the same size again!
    rel_vort.append( np.nan)

    # replace relative vorticity values larger than a limit with the limit value! these results are unphysical and likely due to large jumps
    # in radial distance or total vel.
    # replace with the limit value to still capture the high values in this region
    rel_vort = np.array( rel_vort)

    # only correct values if data are being averaged!
    if window > 1:
        rel_vort[ np.where( rel_vort > .2)[0]] = .2 # np.nan
        # do the same for low values
        rel_vort[ np.where( rel_vort < -.1)[0]] = -.1 # np.nan


    if window > 1:
        # smooth the relative vorticity data to remove peaks!
        rel_vort = pd.Series( rel_vort).rolling(window = window, min_periods=1, center=True).mean()
    # optional code to double smooth the data! this makes it look more like the 2001 paper
    if double_average:
        rel_vort = pd.Series( rel_vort).rolling(window = window, min_periods=1, center=True).mean()

    # correct units (from s-1 to *10^-4 s-1)
    rel_vort = rel_vort * 10000

    return rel_vort


# trim down and create new distance and rmw axes for one flight level dataset!
# inputs:
# yearval: the year of this tc occurence
# fileval: the name of this tc dataset (in the format of "20210929H2_sam.nc")
# tcname: the tc name. helpful for loading track info.
# fl_new: the new flight level dataset. Needed for lat lon values to find distance, and
#         wind speeds / pmins to find RMW wind peaks.
# float_time: the newly created fl time axis (hours, UTC, decimal). Needed for matching
#         track center times to fl times.
# return: the newly created distance and rmw axes, filled with nans if not applicable.
def find_center_dist_rmw( yearval, fileval, tcname, fl_new, float_time):

    #######
    ## step 1.1: see if there's valid track data for this case
    #######
    # make the track file name
    track_data_list = []
    track_file_path = "/Users/etmu9498/research/data/track/"
    trackfilename = tcname + yearval + ".trak"
    # load the files in the track folder: is this tc case present?
    trackfiles = helper_fns.load_flight_level( track_file_path, print_files=False)
    count = trackfiles.count( trackfilename)
    if count == 0:
        print( "No track file for TC " + fileval + " present. Nans saved in dist and rmw axis.")
        nanarray = np.empty( len( float_time))
        nanarray[:] = np.nan
        return nanarray, nanarray
    elif count != 1:
        print( "Duplicate cases?")

    #######
    ## step 1.2: load the tc track data!
    #######
    # account for annoying laura2020 case: extra values given for some reason
    # see if other cases have the same issue??
    track_data_temp = pd.read_fwf( track_file_path + trackfilename, skiprows=3)
    # laura case
    if np.shape( track_data_temp)[1] == 8:
        print('laura error case: update track columns')
        col_names = ['date', 'time (UTC)', 'lat', 'lat dir', 'lon', 'lon dir', 'temp1', 'temp2']
    else:
        col_names = ['date', 'time (UTC)', 'lat', 'lat dir', 'lon', 'lon dir']
    track_data = pd.read_fwf( track_file_path + trackfilename, skiprows=3, names = col_names)

    # add a new column with time as a decimal
    float_times = []
    for ind in range( np.shape( track_data)[0]):
        str_timei = track_data['time (UTC)'][ ind]
        # calculate the decimal time
        h = float( str_timei[0:2])
        m = float( str_timei[3:5])
        s = float( str_timei[6:8])
        float_times.append( h + m / 60 + s / 3600)
    track_data[ 'float_time'] = float_times
    # add a new column with the date in the format of new flight level netcdfs
    newdates = []
    for ind in range( np.shape( track_data)[0]):
        old_datei = track_data.date[ ind]
        newdates.append( old_datei[6:10] + old_datei[0:2]+old_datei[3:5])
    track_data[ 'date2'] = newdates

    #######
    ## step 1.3: trim down the time ranges on the track data to prevent looking at tracks outside the correct limits!
    #######
    # find flight level (fl) time limits for the current and sometimes the next day
    fl_date = fileval[0:8]
    fl_firstt = float_time[0]
    fl_lastt = float_time[-1]
    # cut off all data before the start of the p-3 data
    # find the indices for the times from the first day
    current_day_inds = np.where( track_data['date2'] == fileval[0:8])[0]

    # case where the date of this flight isn't even included in the track list! way outside of bounds.
    # just return nans
    if len( current_day_inds) == 0:
        print( "Flight level date completely before provided track. Nans saved in dist and rmw axis")
        nanarray = np.empty( len( float_time))
        nanarray[:] = np.nan
        return nanarray, nanarray
    # find the closest time to the start time!
    # need to do first_day_inds[0] to account for indices before the start of first_day_inds!
    first_t_ind = current_day_inds[0]  + np.argmin( np.abs( track_data['float_time'][ current_day_inds].values - fl_firstt ))

    newdate = []
    # find the closest time to the end time!
    # same day CASE
    if fl_lastt < 24.0:
        last_t_ind = current_day_inds[0]  + np.argmin( np.abs( track_data['float_time'][ current_day_inds].values - fl_lastt ))
    # new case: the float_time variable shifts into the next day
    # need to look at arrays for the next day
    else:
        # find the next day indices.
        # done this way rather than + 1 to account for flips like 9/29 -> 9/30, or 9/31 -> 10/1
        # look through dates from current day to the end of the array
        for datei in track_data['date2'][ current_day_inds[0] : ]:
            # new date case: save it and end the loop!
            if datei != fileval[0:8]:
                newdate = datei
                break
        # the same sorting if statement as below: find cases that come completely after the last track time :(
        if len( newdate) == 0:
            print( "Flight level date completely after provided track. Nans saved in dist and rmw axis.")
            nanarray = np.empty( len( float_time))
            nanarray[:] = np.nan
            return nanarray, nanarray
        # still in the if statement, find all the inds for this next day
        next_day_inds = np.where( track_data['date2'] == newdate)[0]
        # do - 24 to account for the 24 hour shift!
        last_t_ind = next_day_inds[0]  + np.argmin( np.abs( track_data['float_time'][ next_day_inds].values - ( fl_lastt - 24) ))
    # make the final array!
    dayind = range( first_t_ind, last_t_ind)

    ########
    ## step 1.4: find the tc center distance!
    ## somewhat time consuming but accurate
    ########
    dists = []
    # defining these variables once outside the loop actually saves a ton of time when accessing them
    lat = fl_new.LATref.values
    lon = fl_new.LONref.values
    flnewtime = float_time
    track_lat = track_data.lat.values [dayind]
    track_lon = - track_data.lon.values [dayind] # values are saved as positive from track .txt document
    tracktime = track_data.float_time.values [dayind]
    trackdate = track_data.date2.values [ dayind]

    # case where the date of this flight is included in the track list, but the first time is way past
    # the last time in the track list
    # just return nans
    if len( tracktime) == 0:
        print( "Flight level date completely after provided track. Nans saved in dist and rmw axis.")
        nanarray = np.empty( len( float_time))
        nanarray[:] = np.nan
        return nanarray, nanarray

    # first, account for date changes in the time intervals:
    fl_orig_date = fileval[0:8] # the date when the first pass was completed
    for tracktimei in range( len( tracktime)):
        # new date case: increase the time by 24 hours!
        if trackdate [ tracktimei] != fl_orig_date:
            tracktime[ tracktimei] += 24.0
    latnans = 0
    lonnans = 0
    far_dist = False
    for fl_i in range( len( fl_new.LATref)):
        lati, loni, timei = lat[ fl_i], lon[ fl_i], flnewtime[ fl_i]

        # find the closest center fix time index!
        #########################
        ## 2/13/23 new code: make sure that the closest center fix time is at least smaller than 15 minutes (?)
        ## and that the date is within the correct bounds!
        ## if the time difference is too large, then add nans for this case :/
        ## helpful for avoiding edge case errors ( a 15 minute window is used??)
        #########################
        if np.min( np.abs( tracktime - timei )) > .250 :
            dists.append( np.nan)
            latnans += 1
            lonnans += 1
            # print out a helpful warning to the user for out of range fl times!
            # the far_dist flag is to make sure that this is only printed out once :)
            if not far_dist:
                far_dist = True
                print( "Flight level times are out of range from valid track times (either too early or late.)" +
                        " Nans have been added to the center_dist and rmw array.")
            pass
        else:
            center_i = np.argmin( np.abs( tracktime - timei ))
            center_time = tracktime[ center_i]
            center_date = trackdate[ center_i]
            coords_p3 = np.ma.masked_invalid( ( lati, loni))
            coords_center = np.ma.masked_invalid( ( track_lat[center_i], track_lon[center_i]))
            center_dist_i = distance.geodesic( coords_p3, coords_center).km

            # nan cases- can't find a distance
            if np.isnan( lati):
                latnans += 1
                center_dist_i = np.nan
            if np.isnan( loni):
                lonnans += 1
                center_dist_i = np.nan
            # otherwise, append a valid distance to the dists list
            dists.append( center_dist_i)
        # nice user notices
        if fl_i == 0:
            print( 'number of flight level data points: ' + str( len( fl_new.LATref)))
        if fl_i % 10000 == 0:
            print( 'index = ' + str( fl_i))
    print('number of lat nans: ' + str( latnans))
    print('number of lon nans: ' + str( lonnans))

    #############
    ## next main goal: find an rmw axis!!
    #############
    # step 2.1: find vpeaks using the typical method
    # input vars used later
    max_v_requirement=40
    window=10
    timelim = 60
    spd_avg = pd.Series( fl_new['WS.d']).rolling(window=window, min_periods=1, center=True).mean()
    vpeaks, pmins, time_lims = peak_algs.find_peaks_pmin_time_limit( fl_new['PSURF.d'], spd_avg, float_time, window, timelim=timelim)

    ###########
    # step 2.2:
    ###########
    # find distances at each wind speed peak. Use these distances to define a new rmw axis.
    # go through a couple if statements to determine which rmw1 value to use as the rmw limit
    rmw0_vals = []
    for peak in pmins:
        rmw0_vals.append( dists[ peak])
    rmw1_vals = []
    for peak in vpeaks:
        rmw1_vals.append( dists[ peak])

    # divide up the distance axis into distinct RMW regions
    # do this for every peak
    rmwlist = []
    rainbandlist = []
    for vpeaki, vpeak in enumerate( vpeaks):
        # figure out the distances to the peak winds and current pressure center
        rmw1 = rmw1_vals[ vpeaki]
        # use // to get rid of remainder
        pressurei = vpeaki // 2
        rmw0 = rmw0_vals[ pressurei]
        # first value case: count everything before the first wind speed to the center
        if vpeaki == 0:
            # divide by rmw1 so that the distance at the peak winds = 1!
            rmwlist += (np.array( dists [ 0 : pmins[ pressurei] ]) / rmw1).tolist()
            continue
        # last value case: count everything from the last center to the end of the dataset
        elif vpeaki == len( vpeaks) - 1:
            rmwlist += (np.array( dists [ pmins[ pressurei] : len( dists)]) / rmw1).tolist()
            continue
        # middle cases: things are more complicated lol
        # odd cases: the center has already been passed
        elif vpeaki % 2 == 1:
            # find the index in the rainbands
            # should be halfway between the two closest pressure centers
            # first pressure center + ( difference / 2)
            rainband_ind = int( pmins[ pressurei] + ( ( pmins[ pressurei + 1] - pmins[ pressurei]) / 2) )
            rainbandlist.append( rainband_ind)
            print( 'rainband index = ' + str( rainband_ind))
            # go from the previous pressure center to a defined point in the outer rainbands
            rmwlist += ( np.array( dists [ pmins[ pressurei] : rainband_ind]) / rmw1).tolist()

        # even cases: the center is being approached
        elif vpeaki % 2 == 0:
            # find the index in the rainbands
            rainband_ind = int( pmins[ pressurei - 1] + ( ( pmins[ pressurei ] - pmins[ pressurei - 1]) / 2) )
            rainbandlist.append( rainband_ind)
            print( 'rainband index = ' + str( rainband_ind))
            # go from the past halfway mark out in the rainbands to the next pressure center
            rmwlist += ( np.array( dists [ rainband_ind : pmins[ pressurei] ]) / rmw1).tolist()
        else:
            print( "Error: idk how this happened lol")

    # case with no vpeaks, so not possible to create an rmw axis! just make an empty array for rmw
    if len( rmwlist) == 0:
        nanarray = np.empty( len( float_time))
        nanarray[:] = np.nan
        print( "No wind speed peaks found for this case, so all RMW values will be nans.")
        return dists, nanarray

    # final case: made it through the gauntlet, return the populated dists and rmw arrays!
    return dists, np.array( rmwlist)
