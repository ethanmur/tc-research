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
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023/era5-processing/")
import era_data_processing

# plot type options:
# all_speeds -> plot wind speed and vector directions for both levels
# trimmed_speeds -> plot wind speeds and directions for both levels (4 subplots) within the trimmed tc radius
# 
def wind_plots(year='2021', name = 'henri', lat_spacing = 30, lon_spacing = 30, plot_type='all_speeds', smoothval=1):
    level1 = 250
    level2 = 850
    dpi=150
    cmap='plasma' # jet, plasma
    fs = 12.5

    # load data after getting correct file names for each case
    os.chdir("/Users/etmu9498/research/data/era5/" + year + "/")
    
    # get location and time data for one tc -> for figuring out filenames
    names, dt, all_lats, all_lons = era_data_processing.pull_lat_lons( tc = year)
    tcind = np.where( names == name)[0][0]
    
    # number of cases / era 5 datasets to plot
    ncases = len(dt[tcind])    
    # do this for every case:
    for case in range(ncases):
        datei, lati, loni = dt[tcind][case], all_lats[tcind][case], all_lons[tcind][case]

        latN = lati + lat_spacing # lat needs no changes!
        latS = lati - lat_spacing
        lonW = 360 + loni - lon_spacing    # must be in degrees E (The western hemisphere is captured between 180 and 360 degrees east)
        lonE = 360 + loni + lon_spacing
        # need to reconvert bounds for cartopy plotting
        all_coords = [latN, latS, lonW, lonE]

        bounds = [loni - lon_spacing, loni + lon_spacing, latS, latN]

        # get dates to find file names!
        startdate = datei[0:5] + "/" + year[2:4] + " " + datei[6:8] + ":00" ## must be in this format 'MM/DD/YY HH:MM'
        date = datetime.strptime(startdate,'%m/%d/%y %H:%M')

        # open the correct file
        filestr = 'winds_'+date.strftime('%Y')+date.strftime('%m')+date.strftime('%d')+'_'+date.strftime('%H')+'.nc'
        os.chdir("/Users/etmu9498/research/data/era5/" + year + "/" + name + "/")
        ds = xr.open_dataset(filestr).metpy.parse_cf()

        ######## Read in the variables using a helper function ##########

        # option 1: make simpler wind speed and vector plots for all era 5 data within bounds
        if plot_type == 'all_speeds':
            coords, upperwinds, lowerwinds, spds = era_data_processing.process_era5_data( ds, date, all_coords, find_trimmed_vals=False, smoothval=smoothval)

            ######## Create the upper and lower plots #######
            us = [upperwinds[0], lowerwinds[0]]
            vs = [upperwinds[1], lowerwinds[1]]
            lats, lons = coords[0], coords[1]
            levels = [level1, level2]
            subplot_list = [121, 122]

            #Start the figure and create plot axes with proper projection
            fig = plt.figure(1, figsize=(16, 8))
            helper_fns.change_font_sizes(fs,fs)
            for i in range(2):
                #Set up the projection that will be used
                mapcrs = ccrs.LambertConformal(central_longitude=-55, central_latitude=45, standard_parallels=(33, 45)) 
                mapcrs = ccrs.PlateCarree()

                #Set up the projection of the data; if lat/lon then PlateCarree is what you want
                datacrs = ccrs.PlateCarree()

                ax = plt.subplot(subplot_list[i], projection=mapcrs) 

                print(bounds)
                if len( np.where(np.isnan(bounds))[0]) > 0:
                    print("Non Valid bounds have been provided, likely due to corrupted data. skipping this plot")
                    continue
                    
                ax.set_extent( bounds, ccrs.PlateCarree()) ## Can change the lat/lon bounds

                #Add geopolitical boundaries for map reference
                ax.add_feature(cfeature.LAND, facecolor="#bdbdbd") 
                countries = NaturalEarthFeature(category="cultural", scale="110m", facecolor="none", name="admin_0_boundary_lines_land") 
                ax.add_feature(countries, linewidth=0.5, edgecolor="black") 
                ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5) 
                ax.coastlines('50m', linewidth=0.8)

                #Set up contour and fill intervals
                wspd_levels = arange(0,100,10)
                #Plot the variables
                wspd_fcontours = plt.contourf(lons, lats, spds[i], levels=wspd_levels, cmap=get_cmap('BuPu'), alpha=1, transform=ccrs.PlateCarree())
                # Plot total wind barbs
                udivbrb = us[i] 
                vdivbrb = vs[i]
                wind_slice = (slice(None, None, 10), slice(None, None, 10))
                ax.quiver(lons[wind_slice[0]], lats[wind_slice[1]],
                          udivbrb[wind_slice], # .m
                          vdivbrb[wind_slice], # .m
                          pivot='mid', color='black',
                          #scale=0.5e-6, scale_units='inches',
                          width=0.005,
                          transform=datacrs)

                # add lat / lon labels and ticks!
                gl = ax.gridlines(draw_labels=True) # , linewidth=0)
                gl.top_labels = False
                gl.right_labels = False
                # ax.set_xlabel("") # actually really hard to properly set these labels :(

                #Colorbar and contour labels
                cb = fig.colorbar(wspd_fcontours, orientation='vertical', pad=0.03, extendrect=True, aspect=25, shrink=0.6)
                cb.set_label('Wind Speed (m s$^{-1}$)', size='x-large') 
                plt.grid('off')

                plt.scatter( loni, lati, c='g', s=250, marker='*', label='TC Sam Center', zorder=5)
                leg = plt.legend(loc='upper right', framealpha=1.0, fontsize=fs)
                leg.get_frame().set_linewidth( 1.5) 
                leg.get_frame().set_edgecolor('k')

                plt.title(str( levels[i]) + ' hPa Winds (m s$^{-1}$), ' + startdate, loc='left') 

            print("Figure " + str(case+1) + " created")
                    
            os.chdir("/Users/etmu9498/research/figures/era5-winds/" + year + "/" + name)
            fig.savefig(filestr[: -3] + '_.png', dpi=dpi, bbox_inches='tight')
           
            print("Figure saved.")

        # option 2: make trimmed wind speed and direction plots both levels (4 subplots) for all era 5 data within bounds
        elif plot_type == 'trimmed_speeds':
            tccenter = [lati, loni]
            coords, upperwinds, lowerwinds, spds_trimmed, dirs_trimmed = era_data_processing.process_era5_data( ds, date, all_coords, tccenter, find_trimmed_vals=True, smoothval=smoothval)

            ######## Create the upper and lower plots #######
            us = [upperwinds[0], lowerwinds[0]]
            vs = [upperwinds[1], lowerwinds[1]]
            lats, lons = coords[0], coords[1]
            levels = [level1, level2]
            subplot_list = [221, 222, 223, 224]

            #Set up the projection that will be used
            mapcrs = ccrs.LambertConformal(central_longitude=-55, central_latitude=45, standard_parallels=(33, 45)) 
            mapcrs = ccrs.PlateCarree()
            #Set up the projection of the data; if lat/lon then PlateCarree is what you want
            datacrs = ccrs.PlateCarree()

            #Start the figure and create plot axes with proper projection
            fig = plt.figure(1, figsize=(16, 16))
            for i in range(4):
                ax = plt.subplot(subplot_list[i], projection=mapcrs) 

                if len( np.where(np.isnan(bounds))[0]) > 0:
                    print("Non Valid bounds have been provided, likely due to corrupted data. skipping this plot")
                    continue                    
                ax.set_extent( bounds, ccrs.PlateCarree()) ## Can change the lat/lon bounds

                #Add geopolitical boundaries for map reference
                ax.add_feature(cfeature.LAND, facecolor="#bdbdbd") 
                countries = NaturalEarthFeature(category="cultural", scale="110m", facecolor="none", name="admin_0_boundary_lines_land") 
                ax.add_feature(countries, linewidth=0.5, edgecolor="black") 
                ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5) 
                ax.coastlines('50m', linewidth=0.8)

                #Set up contour and fill intervals
                if np.nanmax( spds_trimmed[0]) > 30:
                    vmin, vmax = 0, 60
                else:
                    vmin, vmax = 0, 30
                # Plot total wind speeds and directions
                if i == 0:
                    spd_label = '250 hPa Winds'
                    wspd_fcontours = plt.pcolormesh(lons, lats, spds_trimmed[0], vmin=vmin, vmax=vmax, cmap=get_cmap(cmap), alpha=1, transform=ccrs.PlateCarree())
                    cb = fig.colorbar(wspd_fcontours, orientation='vertical', pad=0.03, extendrect=True, aspect=25, shrink=0.6)
                    cb.set_label('Wind Speed (m/s)', size='x-large') 
                elif i == 1:
                    wspd_fcontours = plt.pcolormesh(lons, lats, dirs_trimmed[0], vmin=0, vmax=360, cmap=get_cmap(cmap), alpha=1, transform=ccrs.PlateCarree()) 
                    cb = fig.colorbar(wspd_fcontours, orientation='vertical', pad=0.03, extendrect=True, aspect=25, shrink=0.6)
                    cb.set_label('850 hPa Wind Dir. (' + u'\N{DEGREE SIGN}' + ' N)', size='x-large')
                elif i == 2:
                    spd_label = '850 hPa Winds (m/s)'
                    wspd_fcontours = plt.pcolormesh(lons, lats, spds_trimmed[1], vmin=vmin, vmax=vmax, cmap=get_cmap(cmap), alpha=1, transform=ccrs.PlateCarree())
                    cb = fig.colorbar(wspd_fcontours, orientation='vertical', pad=0.03, extendrect=True, aspect=25, shrink=0.6)
                    cb.set_label('Wind Speed (m/s)', size='x-large')                     
                elif i == 3:
                    wspd_fcontours = plt.pcolormesh(lons, lats, dirs_trimmed[1], vmin=0, vmax=360, cmap=get_cmap(cmap), alpha=1, transform=ccrs.PlateCarree())
                    cb = fig.colorbar(wspd_fcontours, orientation='vertical', pad=0.03, extendrect=True, aspect=25, shrink=0.6)
                    cb.set_label('250 hPa Wind Dir. (' + u'\N{DEGREE SIGN}' + ' N)', size='x-large')

                # upper wind arrows
                if i < 2:
                    # Plot total wind barbs
                    udivbrb = upperwinds[0] 
                    vdivbrb = upperwinds[1]
                # lower wind arrows:
                else:
                    udivbrb = lowerwinds[0] 
                    vdivbrb = lowerwinds[1]
                # add wind barbs
                wind_slice = (slice(None, None, 5), slice(None, None, 5))
                ax.quiver(lons[wind_slice[0]], lats[wind_slice[1]],
                          udivbrb[wind_slice], # .m
                          vdivbrb[wind_slice], # .m
                          pivot='mid', color='black',
                          #scale=0.5e-6, scale_units='inches',
                          width=0.005,
                          transform=datacrs)

                # add lat / lon labels and ticks!
                gl = ax.gridlines(draw_labels=True) # , linewidth=0)
                gl.top_labels = False
                gl.right_labels = False
                gl.xlines = False
                gl.ylines = False
                # ax.set_xlabel("") # actually really hard to properly set these labels :(


                # add tc center location
                plt.scatter( loni, lati, c='g', s=250, marker='*', label='TC ' + name.title() + ' Center', zorder=5)
                # add a fake arrow for the legend!
                testarrow = plt.scatter( -1000, -1000, c='k', marker=r'$\longrightarrow$', s=150, label=spd_label )
                leg = plt.legend(loc='upper right', framealpha=1.0, fontsize=fs)
                leg.get_frame().set_linewidth( 1.5) 
                leg.get_frame().set_edgecolor('k')

                plt.title(spd_label + ', ' + startdate, loc='left') 

            print("Figure " + str(case+1) + " created\n")
                    
            os.chdir("/Users/etmu9498/research/figures/era5-winds-trimmed/" + year + "/" + name)
            fig.savefig(filestr[: -3] + '_.png', dpi=dpi, bbox_inches='tight')



def shear_comparison_plot(year='2021', name = 'henri', manual_name=False):

    if year != '2021':
        print("Please only use year=2021 for this analysis!")
        return

    # get ships data via helper function here
    dt, shearmag, sheardir = era_data_processing.pull_shear( name=name)
    # convert shearmag from kt to m/s
    shearmag = np.array( shearmag) * 0.514444

    print(len(dt))

    print(dt[1])
    print(type(dt[1]))

    # get era5 derived shear here! saved locally after a fairly long calculation
    os.chdir("/Users/etmu9498/research/data/era5-shear/")
    # load a special dataset (smoothed, etc) with this flag
    if manual_name:
        nc_shear = xr.open_dataset( manual_name)
    else:
        nc_shear = xr.open_dataset( "shear_" + name + ".nc")


    # datenew = []
    # datestrings = np.datetime_as_string(nc_shear.date.values)
    # for i, val in enumerate( datestrings):
    #     datenew.append( val[5:7]+ '/' + val[8:10] + ' ' + )
    datenew = nc_shear.date.values

    print(datenew[1])
    print(type(datenew[1]))

    plt.figure( figsize = (12, 8))
    helper_fns.change_font_sizes( 14, 14)
    lw = 2
    n = 2  # Keep every 4th label

    plt.subplot( 211)
    plt.title( "SHIPS Derived Plots for TC " + name.title() + ", " + year)
    
    plt.plot( dt, shearmag, c='r', linewidth=lw, label='SHIPS Data')
    plt.plot( datenew, nc_shear.shear_spd, c='b', label='ERA 5 Data')

    if manual_name:
        plt.plot( datenew, nc_shear.shear_avg, c='y', label='ERA 5 Average')


    plt.ylabel( "Shear Mag (m/s)")
    ax = plt.gca()
    [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
    # plt.xticks(rotation=45, ha="right")
    plt.grid()
    plt.legend(loc='upper right')

    plt.subplot( 212)

    # correct simple delta shear direction here! If a negative value is found, just wrap it around to 360 degrees again?
    shear_change = nc_shear.shear_dir_change.values
    negative_inds = np.where( shear_change < 0.0)[0]
    shear_change[negative_inds] = 360 + shear_change[negative_inds]

    plt.plot( dt, sheardir, c='r', linewidth=lw, label='SHIPS Data')
    plt.plot( datenew, nc_shear.shear_dir_avg, c='b', linewidth=lw, label='ERA 5 Avg')
    # plt.plot( datenew, shear_change, c='y', linewidth=lw, label='ERA 5 Change')
    
    plt.ylabel( "Shear Dir (Degrees N)")
    # plt.xticks(rotation=45, ha="right")
    plt.grid()
    ax = plt.gca()
    [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
    # plt.xlabel( "Hours Since " + starttime + " UTC on " + startdate)
    plt.xlabel( "Date")

    fig=plt.gcf()
    fig.tight_layout()
    plt.legend(loc='upper right')
