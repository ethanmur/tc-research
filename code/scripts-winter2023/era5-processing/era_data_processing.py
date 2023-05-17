import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import metpy.calc as mpcalc
from metpy.units import units
from numpy import *
import xarray as xr
import pandas as pd
import numpy as np
from netCDF4 import Dataset, num2date
import math
import sys
from datetime import datetime, timedelta

import cdsapi
import os
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns



# using a combination of ships data (for 2021 cases) and nhc tc reports (for 2022 cases),
# figure out which dates have active TC events, and get the lat / lon tc centers for these dates! where the shear centroid
# will be centered
def pull_lat_lons( tc = '2021'):
    
    # process data from SHIPS dataset to set up for 2021 intensification analysis
    # this code is based on "/code/eye-cloud-paper/figure 4 cloud dists vs TC intensification.ipynb"
    if tc == '2021':
        # define tc names used in this analysis
        tcnames_orig = ['fred', 'grace', 'henri', 'ida', 'sam']
        tcnames = []
        # only keep first 4 letters to match the format of the ships data
        for i in range(len( tcnames_orig)):
            tcnames.append( tcnames_orig[i][0:4].upper())
        year = '2021'

        # load the full ships dataset
        print( "Loading SHIPS dataset")
        os.chdir(  "/Users/etmu9498/research/data/ships/")
        file1 = open('lsdiaga_1982_2021_sat_ts_5day.txt', 'r')
        Lines = file1.readlines()
        # there are a bunch of lines in the total ships file!!
        print( "SHIPS dataset created")
        print( "Number of lines in dataset: " + str( len( Lines)))

        # lists of lists holding important date, intensity, etc information for each tc of interest!
        all_dates = []
        all_time_since_start = []
        all_datetimes = []
        all_lat = []
        all_lon = []
        all_startdate = []
        all_starttime = []

        # do this for each tc name listed above
        for strname in tcnames:
            # convenient variables for plotting
            dates = []
            time_since_start = []
            datetimes = []
            lat = []
            lon = []
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
                        datetime_orig = datetime( int( year), month, day, hours)
                        datetimes.append( datetime_orig.strftime( "%m/%d %Hh"))
                        break
                # search for lat and lon values!
                for i in range( headerval,  len( Lines) ):
                    if 'LON' in Lines[i]:
                        lon.append( - float( Lines[i][12 : 15]) / 10 ) # the last 3 vals for vmax at 0 hours
                        break
                for i in range( headerval,  len( Lines) ):
                    if 'LAT' in Lines[i]:
                        lat.append( float( Lines[i][11:15] ) / 10.  ) # the last 3 vals for vmax at 0 hours
                        break

            # once done looping through the headers, append the list of dates, pressures, etc to the all_ lists
            all_dates.append( dates)
            all_time_since_start.append( time_since_start)
            all_datetimes.append( datetimes)
            all_lat.append( lat)
            all_lon.append( lon)
            all_startdate.append( startdate)
            all_starttime .append( starttime)
        return np.array(tcnames_orig), all_datetimes, all_lat, all_lon
            
    # other case: process munged data from NHC reports to set up for 2022 intensification analysis
    elif tc == '2022':
        tcnames_orig_2022 = ['earl', 'fiona', 'ian', 'julia']
        # months aren't saved in the original excel files :/ so add them here! from noaa summary docs
        months = ['09', '09', '09', '10']

        # lists of lists holding important date, intensity, etc information for each tc of interest!
        all_datetimes = []
        all_lons = []
        all_lats = []

        # where 2022 munged data resides
        data_path_2022 = "/Users/etmu9498/research/data/intensity-6-hour-updates-temp/"

        # do this for each tc name listed above
        for stri, strname in enumerate( tcnames_orig_2022):
            datetimes = []

            # read in this tc's dataset
            os.chdir(data_path_2022)
            data = pd.read_excel( "tc-" + strname + ".xlsx")

            # cycle through each datetime and correct it
            for di, dval in enumerate(data['Date/Time (UTC)'].values):
                datetimes.append( months[stri] + '/' + dval[0:2] + ' ' + dval[5:7] + 'h')

            # append the list of dates, pressures, etc to the all_ lists
            all_datetimes.append( datetimes )
            all_lats.append( data['Latitude (ON)'].values.tolist())
            all_lons.append( (-1 * data['Longitude (OW)'].values).tolist())
            
        return np.array(tcnames_orig_2022), all_datetimes, all_lats, all_lons



# download netcdf era 5 data for a particular tc here! do the analysis in a later step
# enter a suitable year and tc name to download data
def download_era5( year='2021', name = 'henri'):

    ############# Declare time, level, and lat/lon boundaries here ##########
    # do so automatically! will be useful code later :)
    # inputs
    level1 = 250
    level2 = 850

    # get location and time data for one tc
    names, dt, all_lats, all_lons = pull_lat_lons( tc = year)
    tcind = np.where( names == name)[0][0]
    
    # number of cases / era 5 datasets to download
    ncases = len(dt[tcind])
    
    # do this for every case:
    for case in range(ncases):
        datei, lati, loni = dt[tcind][case], all_lats[tcind][case], all_lons[tcind][case]

        # convert saved dates, etc into format for era 5 analysis!
        startdate = datei[0:5] + "/" + year[2:4] + " " + datei[6:8] + ":00" ## must be in this format 'MM/DD/YY HH:MM'
        date = datetime.strptime(startdate,'%m/%d/%y %H:%M')

        ######### Download ERA5 data using CDS-API ##########
        filestr = 'winds_'+date.strftime('%Y')+date.strftime('%m')+date.strftime('%d')+'_'+date.strftime('%H')+'.nc'     # you can change these filenames as you see fit

        os.chdir("/Users/etmu9498/research/data/era5/" + year + "/" + name + "/")
        c = cdsapi.Client()
        c.retrieve('reanalysis-era5-pressure-levels', {
            'product_type': 'reanalysis',
            'variable': ['u_component_of_wind', 'v_component_of_wind'],
            'pressure_level': ['250', '850'],
            'year': date.strftime('%Y'),         # substitutes the year based off the startdate in the above cell
            'month': date.strftime('%m'),        # substitutes the month based off the startdate in the above cell
            'day': date.strftime('%d'),          # substitutes the year based off the startdate in the above cell
            'time': date.strftime('%H:%M'),      # substitutes the time based off the startdate in the above above
            'format': 'netcdf',
            }, filestr)
           
        print(filestr)
        print("Dataset " + str(case+1) + "/" + str(ncases) + ' downloaded')


# inputs: ds, a properly loaded era5 dataset in .netcdf format
def process_era5_data( ds, date, all_coords, smoothval=1, find_trimmed_vals=False):
        level1, level2 = 250, 850

        latN, latS, lonW, lonE = all_coords[0], all_coords[1], all_coords[2], all_coords[3]

        lats = ds.latitude.metpy.sel(latitude=slice(latN,latS))
        lons = ds.longitude.metpy.sel(longitude=slice(lonW,lonE))
        # upper wind
        uwnd1 = ds.u.metpy.sel(time=date,level = level1,latitude=slice(latN,latS),longitude=slice(lonW,lonE))
        vwnd1 = ds.v.metpy.sel(time=date,level = level1,latitude=slice(latN,latS),longitude=slice(lonW,lonE))
        # lower wind
        uwnd2 = ds.u.metpy.sel(time=date,level = level2,latitude=slice(latN,latS),longitude=slice(lonW,lonE))
        vwnd2 = ds.v.metpy.sel(time=date,level = level2,latitude=slice(latN,latS),longitude=slice(lonW,lonE))
        ######### Optional: Smooth the variables ##########
        uwnd1 = mpcalc.smooth_gaussian(uwnd1,smoothval)
        vwnd1 = mpcalc.smooth_gaussian(vwnd1,smoothval)
        uwnd2 = mpcalc.smooth_gaussian(uwnd2,smoothval)
        vwnd2 = mpcalc.smooth_gaussian(vwnd2,smoothval)

        wspd1 = mpcalc.wind_speed(uwnd1,vwnd1)   # wind speed
        wspd2 = mpcalc.wind_speed(uwnd2,vwnd2)

        # default case: just return the positions and global wind speeds... don't trim by radial distance!
        if not find_trimmed_vals:
            return (lats, lons), (uwnd1, vwnd1), (uwnd2, vwnd2), (wspd1, wspd2)

        # extra processing case: trim wind speeds and directions down to the 200 to 800 km circle
        else:
            # create a distance from the TC center axis!
            # save distances here
            center_dists = np.empty( np.shape(uwnd1))
            for i in range(np.shape(uwnd1)[0]):
                lat = lats[i]
                for j in range(np.shape(uwnd1)[0]):
                    lon = lons[j]
                    center_dists[i][j] = distance.geodesic( center, [lat, lon]).km
                if i % 20 == 0:
                    print("Case " + str(i) + " complete")        
                    
            # replace distances with nans if it's smaller than 200 km, larger than 800 km
            center_dists_donut = np.empty( np.shape(uwnd1))
            # search down columns
            for i in range(np.shape(center_dists_donut)[0]):
                # search through rows
                for j in range(np.shape(center_dists_donut)[1]):
                    if center_dists[i][j] >= 200 and center_dists[i][j] <= 800:
                        center_dists_donut[i][j] = center_dists[i][j]
                    else:
                        center_dists_donut[i][j] = np.nan
                # if i % 20 == 0:
                #     print("Case " + str(i) + " complete")

                # find wind speed and direction at each grid point using metpy
                wspd1 = mpcalc.wind_speed(uwnd1,vwnd1)   # wind speed
                wspd2 = mpcalc.wind_speed(uwnd2,vwnd2)
                wdir1 = mpcalc.wind_direction(uwnd1,vwnd1, convention='from')
                wdir2 = mpcalc.wind_direction(uwnd2,vwnd2, convention='from')

                # replace wind directions and speeds with nans if it's smaller than 200 km, larger than 800 km from tc center
                wspd1_trim = np.empty( np.shape(uwnd1))
                wspd2_trim = np.empty( np.shape(uwnd1))
                wdir1_trim = np.empty( np.shape(uwnd1))
                wdir2_trim = np.empty( np.shape(uwnd1))
                # search down columns
                for i in range(np.shape(center_dists_donut)[0]):
                    # search through rows
                    for j in range(np.shape(center_dists_donut)[1]):
                        if center_dists[i][j] >= 200 and center_dists[i][j] <= 800:
                            wspd1_trim[i][j] = wspd1[i][j]
                            wspd2_trim[i][j] = wspd2[i][j]
                            wdir1_trim[i][j] = wdir1[i][j]
                            wdir2_trim[i][j] = wdir2[i][j]
                        else:
                            wspd1_trim[i][j] = np.nan
                            wspd2_trim[i][j] = np.nan
                            wdir1_trim[i][j] = np.nan
                            wdir2_trim[i][j] = np.nan

            return (lats, lons), (uwnd1, vwnd1), (uwnd2, vwnd2), (wspd1_trim, wspd2_trim), (wdir1_trim, wdir2_trim)