# import...
import numpy as np
import os
import sys
import pandas as pd
from geopy import distance
import metpy.calc as mpcalc
import xarray as xr
import math
import matplotlib.pyplot as plt
os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops

# code below taken from "tests/2023-06-10 shear quadrant height dists",
# "scripts-winter2023/cloud-top-height-stats/distributions_general.py", or 
# "code/tests/2023-06-07 find cartesian and shear xy distances.ipynb"

# Goal: create a generalized function that finds cloud heights and compiles info for all composite groups (intensity, shear, etc).
#       do all the sorting later: the information should all be in this dataframe!
#       can also save locally! easier access later
#       code is partially based on "code/eye cloud paper figures/Figure 4 cloud dists vs shear quadrant.ipynb"
# Inputs: None! Just run this script to create the nice, organized data file
#        No need for a tc='all' input: just do the sorting after making this dataframe!
# Return: a dataframe with height / distribution information for every pass! 
# Notes * -> a list, category = WH, intensifying, DSR depending on input flags! Structure:
# flight | pass | UTC time * | x dists * | y dists * | cloud heights * | vertical bins | normalized dists | defined eyewalls or not | intensity | intensification | tc category | shear strength | shear dir
# case 0
# case 1
# ... 
def find_dists(save=True):
    # use a helper fn to get year and file names, along with intensity, etc. metadata
    tc='all'
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc=tc)
    metadata = eyewall_metadata.all_metadata()
    
    # the dataframe saved below contains simple xy distances for each day. 
    df_dists = shear_date_setup()
    # use the helper function below to go from single dates to individual eye passes. Add to this dataframe for the rest of the code
    df_dists = separate_passes( df_dists)
    # add intensity, shear, category, etc info during this step!
    df_final = add_info( df_dists)

    if save:
        savepath = "/Users/etmu9498/research/data/aa_paper_1_data"
        df_dists.to_pickle(savepath + "/flight_heights_metadata.pkl")
    return df_final



# use this step to create the initial dataframe holding date, pass, and position info for every
def shear_date_setup():
    # use a helper fn to get the relevant years and files
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc='all')

    alldates = []
    allxdists = []
    allydists = []
    alltimes = []
    
    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( crlfilelist)):
        # count all the names in this year, and add to the count
        filecount += len( crlfilelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))
    print( 'data saved to: /Users/etmu9498/research/figures/CRL-all-data-processed/"year"-clusters-"cluster type')

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):
            print("Date: " + crlname[3:13])
                        
            os.chdir("/Users/etmu9498/research/data/crl-all-data-processed/" + yearval + "/")
            crl_data = xr.open_dataset(crlname)
            float_time = crl_data.time.values
            
            # print(float_time[0])
            # print(float_time[-1])
            
            # step 0: 
            # find a matching flight level dataset! Given the year and crl file name
            # first, go to the correct year folder for the crl / fl case
            flpath = "/Users/etmu9498/research/data/in-situ-noaa-processed/" + yearval
            # get a list of the files within this folder
            fllist = helper_fns_winter2023.load_flight_level( flpath, print_files=False)
            # go through all of the files within fllist. Save all the ones with matching dates
            # with the crl data in a new list!
            matchlist = []
            for fl_filei, fl_fileval in enumerate( fllist):
                if crlname[3:-13] == fl_fileval[ : 10]:
                    matchlist.append( fl_fileval)
            # error cases: either too many files with the same date from the same plane, or no valid files!
            if len( matchlist) == 0:
                print(crlname[3:-18])
                print( "Error! No flight level dataset found. Investigate!")
            elif len( matchlist) > 1:
                print( "Error! Too many flight level datasets found for the same date and aircraft. Investigate!")
                print(crlname[3:-18])
            else:
                flname = matchlist[0]
                tcname = flname[ 11 : -13]
            
            #######
            ## step 1.1: see if there's valid track data for this case
            #######
            # make the track file name
            track_data_list = []
            track_file_path = "/Users/etmu9498/research/data/track/"
            trackfilename = tcname + yearval + ".trak"
            # load the files in the track folder: is this tc case present?
            trackfiles = helper_fns_winter2023.load_flight_level( track_file_path, print_files=False)
            count = trackfiles.count( trackfilename)
            if count == 0:
                print( "No track file for TC " + tcname + " present. Nans saved in dist and rmw axis.")
                nanarray = np.empty( len( float_time))
                nanarray[:] = np.nan
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
            fl_date = flname[0:8]
            fl_firstt = float_time[0]
            fl_lastt = float_time[-1]
            # cut off all data before the start of the p-3 data
            # find the indices for the times from the first day
            current_day_inds = np.where( track_data['date2'] == flname[0:8])[0]

            # case where the date of this flight isn't even included in the track list! way outside of bounds.
            # just return nans
            if len( current_day_inds) == 0:
                print( "Flight level date completely before provided track. Skipping this case.")
                nanarray = np.empty( len( float_time))
                nanarray[:] = np.nan

                continue
                
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
                    if datei != flname[0:8]:
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
            lat = crl_data.Lat.values
            lon = crl_data.Lon.values
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
                continue
         
            # first, account for date changes in the time intervals:
            fl_orig_date = flname[0:8] # the date when the first pass was completed
            for tracktimei in range( len( tracktime)):
                # new date case: increase the time by 24 hours!
                if trackdate [ tracktimei] != fl_orig_date:
                    tracktime[ tracktimei] += 24.0
            latnans = 0
            lonnans = 0
            far_dist = False
            
            xdists, ydists = [], []
            for fl_i in range( len( crl_data.Lat)):
                lati, loni, timei = lat[ fl_i], lon[ fl_i], float_time[ fl_i]

                # find the closest center fix time index!
                #########################
                ## 2/13/23 new code: make sure that the closest center fix time is at least smaller than 15 minutes (?)
                ## and that the date is within the correct bounds!
                ## if the time difference is too large, then add nans for this case :/
                ## helpful for avoiding edge case errors ( a 15 minute window is used??)
                #########################
                if np.min( np.abs( tracktime - timei )) > .50 :
                    latnans += 1
                    lonnans += 1
                    # print out a helpful warning to the user for out of range fl times!
                    # the far_dist flag is to make sure that this is only printed out once :)
                    if not far_dist:
                        far_dist = True
                        print( "Flight level times are out of range from valid track times (either too early or late.)" +
                                " Nans have been added to the center_dist and rmw array.")
                    continue
                else:
                    center_i = np.argmin( np.abs( tracktime - timei ))
                    center_time = tracktime[ center_i]
                    center_date = trackdate[ center_i]
                    coords_p3 = np.ma.masked_invalid( ( lati, loni))
                    coords_center = np.ma.masked_invalid( ( track_lat[center_i], track_lon[center_i]))

                    # center_dist_i = distance.geodesic( coords_p3, coords_center).km                    
                    # xdist_simplei =  ( abs(coords_center[1]) - abs(coords_p3[1])) * 111.32 * math.cos( math.radians(coords_p3[0]))
                    # ydist_simplei =  (coords_p3[0] - coords_center[0]) * 110.574
                    
                    xdist_simplei = distance.geodesic( np.ma.masked_invalid((track_lat[center_i], loni)), coords_center).km # .astype(str).str[:-3].astype(float)
                    ydist_simplei = distance.geodesic( np.ma.masked_invalid((lati, track_lon[center_i])), coords_center).km # .astype(str).str[:-3].astype(float)
                    
                    # - x correction case
                    if abs( loni) > abs( track_lon[center_i]): 
                        xdist_simplei = - xdist_simplei
                    # - y case
                    if abs( lati) < abs( track_lat[center_i]):
                        ydist_simplei = - ydist_simplei
                    
                    # nan cases- can't find a distance
                    if np.isnan( lati):
                        latnans += 1
                        xdist_simplei = np.nan
                        ydist_simplei = np.nan
                    if np.isnan( loni):
                        lonnans += 1
                        xdist_simplei = np.nan
                        ydist_simplei = np.nan
                    # otherwise, append a valid distance to the dists list
                    xdists.append( xdist_simplei)
                    ydists.append( ydist_simplei)
      
            # print out helpful notices and add distances to total list!
            print('number of lat nans: ' + str( latnans))
            print('number of lon nans: ' + str( lonnans))
            allxdists.append(xdists)
            allydists.append(ydists)
            alltimes.append(float_time)
            alldates.append(crlname[3:13])

    df_dists = pd.DataFrame()
    df_dists['dates'] = alldates
    df_dists['times'] = alltimes
    df_dists['xdists'] = allxdists
    df_dists['ydists'] = allydists
    return df_dists



# Helper fn used to take xy distances from individual days and split them into specific passes!
# inputs: df_dists, a pandas dataframe with date, time, and xy distance info for single days.
# return: df_new, a pandas dataframe with original and shear corrected xy distances for each individual eye pass.
# single layer cloud data are used for this analysis: easier to match single heights to single xy distance
def separate_passes( df_dists):
    # the same code as above, but only do this for TC eye passes!
    tc = 'all'
    multi_layers = False
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc='all')
    metadata = eyewall_metadata.all_metadata()
    crl_root_path = "/Users/etmu9498/research/data/crl-all-data-processed/"
    # add trimmed time and distance values here
    datetrim, passtrim, timetrim, xtrim, ytrim, xsheartrim, ysheartrim, cloudheights, definedlist = [], [], [], [], [], [], [], [], []
    # 8/7/23 new step: append rmw and radial distances- should already be saved in the processed crl datasets
    radialdists = []
    rmws = []

    alltimes = df_dists['times'].values
    allxdists = df_dists['xdists'].values
    allydists = df_dists['ydists'].values
    shearxdists, shearydists = find_shear_dists(df_dists) # use helper fn to find shear shifted xy distances

    # keep track of which tc we're on with this counter
    flightcounter = 0
    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):    
            # get eye time limits for each date
            date = crlname[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if crlname[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                    defined_eyewalls = metadata[ yearval]['defined_eyewall'][ '0812am']
                elif crlname[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
                    defined_eyewalls = metadata[ yearval]['defined_eyewall'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
                defined_eyewalls = metadata[ yearval]['defined_eyewall'][ date]
            else:
                eyewall_limits = [ ()]
            # do this for each of the eyewall limit pairs! Can have multiple eyes per crl dataset
            for eyei, eyeval in enumerate( eyewall_limits):
                if len( eyeval) > 0:
                    # find the corresponding indices to the time limits
                    ind0 = np.argmin( np.abs( alltimes[flightcounter] - eyeval[0] ))
                    ind1 = np.argmin( np.abs( alltimes[flightcounter] - eyeval[1] ))
                    # clip relevant fields down to the eyewall limits
                    # load crl data
                    os.chdir( crl_root_path + yearval)
                    crl_data = xr.open_dataset( crlname)
                    H = crl_data.height
                    power = crl_data.P_ch1[ ind0 : ind1, :]
                    axis = crl_data.time[ ind0 : ind1]
                    p3_height = crl_data.p3_height[ ind0 : ind1]
                    # find cloud top heights for values within the specified eye distance range
                    if yearval == '2021':
                        cutoff = -30
                    elif yearval == '2022':
                        cutoff = -40
                    if multi_layers:
                        heights_lists, time, count = find_cloud_tops.find_multi_cloud_heights( H, power, axis, p3_height, cutoff_power = cutoff)
                        heights = np.array( [item for sublist in heights_lists for item in sublist])
                    else:
                        # regular case: find the top cloud height layer
                        heights, time = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = cutoff)
                    # add xy and height data to the proper lists                    
                    if ind0 != ind1:
                        datetrim.append( crlname[3:13])
                        passtrim.append( eyei)
                        xtrim.append( allxdists[flightcounter][ind0:ind1])
                        ytrim.append( allydists[flightcounter][ind0:ind1])
                        xsheartrim.append( shearxdists[flightcounter][ind0:ind1])
                        ysheartrim.append( shearydists[flightcounter][ind0:ind1])
                        cloudheights.append( heights)
                        timetrim.append( alltimes[flightcounter][ind0:ind1])
                        definedlist.append( defined_eyewalls)    
                        radialdists.append(crl_data.center_dist.values[ind0:ind1])
                        rmws.append(crl_data.rmw.values[ind0:ind1])

                    else:
                        print("Error: matching indices :(")
            # skip over two dates that create problems for some reason
            if crlname == 'P3_20220831H1_processed.nc' or crlname == 'P3_20220830H1_processed.nc':
                pass # print('ERROR CASE!!')
            else:
                flightcounter += 1
    # convert lists to a nice dataframe!
    df_dists_trim = pd.DataFrame()
    df_dists_trim['flight'] = datetrim
    df_dists_trim['pass'] = passtrim
    df_dists_trim['times'] = timetrim
    df_dists_trim['radial dists'] = radialdists
    df_dists_trim['rmw'] = rmws
    df_dists_trim['xdists'] = xtrim
    df_dists_trim['ydists'] = ytrim
    df_dists_trim['xdistsshear'] = xsheartrim
    df_dists_trim['ydistsshear'] = ysheartrim
    df_dists_trim['cloudheights'] = cloudheights
    df_dists_trim['Defined Eyewalls'] = definedlist
    return df_dists_trim


# after creating a dataframe with passes, cloud heights, and distances in separate columns, 
# add more info like shear and intensity here. Similar methods as before!
def add_info( df_dists):
    # the same code as above, but only do this for TC eye passes!
    tc = 'all'
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc=tc)
    metadata = eyewall_metadata.all_metadata()
    intensity_l, intensification_l, category_l, shearmag_l, sheardir_l = [], [], [], [], []

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):        
            # get eye time limits for each date
            date = crlname[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if crlname[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                    intensity = metadata[ yearval]['intensity'][ '0812am']
                    intensification = metadata[ yearval]['intensification'][ '0812am']
                    category = metadata[ yearval]['category'][ '0812am']
                    shearmag = metadata[ yearval]['shear_mag'][ '0812am']
                    sheardir = metadata[ yearval]['shear_dir'][ '0812am']
                elif crlname[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
                    intensity = metadata[ yearval]['intensity'][ '0812pm']
                    intensification = metadata[ yearval]['intensification'][ '0812pm']
                    category = metadata[ yearval]['category'][ '0812pm']
                    shearmag = metadata[ yearval]['shear_mag'][ '0812pm']
                    sheardir = metadata[ yearval]['shear_dir'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
                intensity = metadata[ yearval]['intensity'][ date]
                intensification = metadata[ yearval]['intensification'][ date]
                category = metadata[ yearval]['category'][ date]
                shearmag = metadata[ yearval]['shear_mag'][ date]
                sheardir = metadata[ yearval]['shear_dir'][ date]
            else:
                eyewall_limits = [ ()]        
            # do this for each of the eyewall limit pairs! Can have multiple eyes per crl dataset
            for eyei, eyeval in enumerate( eyewall_limits):
                if ~ np.isnan(intensity):
                    intensity_l.append( intensity)
                    intensification_l.append( intensification)
                    category_l.append( category)
                    shearmag_l.append( shearmag)
                    sheardir_l.append( sheardir)
    # add the new categories to the dataframe!
    df_dists['intensity'] = intensity_l
    df_dists['intensification'] = intensification_l
    df_dists['category'] = category_l
    df_dists['shearmag'] = shearmag_l
    df_dists['sheardir'] = sheardir_l
    return df_dists


# helper fns (x and y are arrays)
def cart2pol(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return(r, theta)
def pol2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return(x, y)

# go from cartesian xy to shear shifted xy distances here
def find_shear_dists( df_dists):
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc='all')
    metadata = eyewall_metadata.all_metadata()
    allxdists = df_dists['xdists'].values
    allydists = df_dists['ydists'].values

    sheardirlist = []
    # pulling shear for each case- make sure not to get shear from 08/30 and 08/31 which have no data!
    # keep track of which tc we're on with this counter
    flightcounter = 0
    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):        
            # get eye time limits for each date
            date = crlname[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if crlname[11:13] == "H1":
                    sheardir = metadata[ yearval]['shear_dir'][ '0812am']
                elif crlname[11:13] == "H2":
                    sheardir = metadata[ yearval]['shear_dir'][ '0812pm']
            elif date in metadata[ yearval]['shear_dir'].keys():
                sheardir = metadata[ yearval]['shear_dir'][ date]
            else:
                sheardir = [ ()]            
            # deal with problem cases
            if crlname == 'P3_20220831H1_processed.nc' or crlname == 'P3_20220830H1_processed.nc':
                pass # print('ERROR CASE!!')
            else:
                sheardirlist.append( sheardir)
                flightcounter += 1
    # save new x and y dists below
    shearxdists, shearydists = [], []
    # do this for every set of x and y value
    for i in range(len(sheardirlist)):
        r, theta = cart2pol( np.array(allxdists[i]), np.array(allydists[i]))        
        theta2 = theta + np.radians( sheardirlist[i])       
        xtemp, ytemp = pol2cart(r, theta2)
        shearxdists.append(xtemp)
        shearydists.append(ytemp)
    return shearxdists, shearydists
