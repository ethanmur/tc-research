"""
Created 6/1/22
Last edited 6/2/22

Description: this file contains scripts used to automatically generate the
             images used to create animated gifs of tc motion, using a
             variety of GOES satellite bands.

"""

# import packages used by these scripts
from datetime import datetime
import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy  # noqa: F401
import numpy as np
import xarray as xr
import shapely.geometry as sgeom
import matplotlib.patches as mpatches
from subprocess import run

os.chdir("/Users/etmu9498/research/code/scripts")
import helper_fns

def load_goes( gif_path):
    return helper_fns.display_data_files( gif_path, 'goes', 'hide-list')


def goes_ir_gif( file_names, gif_path, crl_name, crl_path): # flight_line):
    """
    goes_ir_gif generates images, whose data is located in a folder specified
    by the user, and saves them automatically to a unique folder. The output
    folder can be found under research/figures/goes-gifs

    :param
    :param flight_line: Should this function plot the P-3 flight path above the
                        IR image? Inputs can be either true or false.
    :return: None. The function implicitly saves images to an automatically
             generated folder.
    """

    # make sure that there's a folder in place to hold output images!
    os.chdir( gif_path)
    # use the date from the first image to create a new folder
    first_dataset = xr.open_dataset( file_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    # check if folder already exists. If not, create folder
    os.chdir( "/Users/etmu9498/research/figures/goes-gifs/")
    if not os.path.isdir( str( scan_start_first_dataset.strftime('%m%d%Y')) + "_ir"):
        os.makedirs( "/Users/etmu9498/research/figures/goes-gifs/" + str( scan_start_first_dataset.strftime('%m%d%Y')) + "_ir")
        print( 'New folder created')
    else:
        print( 'Existing folder accessed')
    os.chdir(gif_path)

    # figure out the maximum and minimum brightness temperatures in all datasets
    # to make the colorbar constant between runs
    # just some initial valus that will be automatically rewritten
    t_max = 0
    t_min = 1000
    for goes_ind in range( len( file_names)):
        C = xr.open_dataset( file_names[ goes_ind])
        if np.max( C['CMI_C13'].data) > t_max:
            t_max = np.max( C['CMI_C13'].data )
        if np.min( C['CMI_C13'].data) < t_min:
            t_min = np.min( C['CMI_C13'].data )

    # make images for gif
    for goes_ind in range( len( file_names)):
        # load data for this specific goes_ind
        C = xr.open_dataset( file_names[ goes_ind])

        # Scan's start time, converted to datetime object
        # Scan's end time, converted to datetime object
        # File creation time, convert to datetime object
        # The 't' variable is the scan's midpoint time
        scan_start = datetime.strptime(C.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
        scan_end = datetime.strptime(C.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
        file_created = datetime.strptime(C.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')
        midpoint = str(C['t'].data)[:-8]
        scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

        """
        # applying a normalization filter to channel 13 to correct IR stuff!
        cleanIR = C['CMI_C13'].data
        # Normalize the channel between a range.
        cleanIR = (cleanIR-90)/(313-90)
        # Apply range limits to make sure values are between 0 and 1
        cleanIR = np.clip(cleanIR, 0, 1)
        # Invert colors so that cold clouds are white
        cleanIR = 1 - cleanIR
        # Lessen the brightness of the coldest clouds so they don't appear so bright
        # when we overlay it on the true color image.
        cleanIR = cleanIR/1.4
        stacked_cleanIR = np.dstack([cleanIR, cleanIR, cleanIR])
        """

        print("Starting to make image " + str( goes_ind + 1))

        # make figure
        fig = plt.figure(figsize=(15, 12))
        # these steps use goes ch 2 as a proxy to load things like lat and lon positions
        dat = C.metpy.parse_cf('CMI_C02')
        x = dat.x
        y = dat.y
        geos = dat.metpy.cartopy_crs
        ax = fig.add_subplot(1, 1, 1, projection=geos)

        """
        max_lon = C.geospatial_lat_lon_extent.geospatial_westbound_longitude
        min_lon = C.geospatial_lat_lon_extent.geospatial_eastbound_longitude
        max_lat = C.geospatial_lat_lon_extent.geospatial_northbound_latitude
        min_lat = C.geospatial_lat_lon_extent.geospatial_southbound_latitude
        """

        # plot goes data
        # 2.5 * stacked_cleanIR
        img = ax.imshow( C['CMI_C13'].data, origin='upper', cmap='Greys',
                  extent=(x.min(), x.max(), y.min(), y.max()),
                  transform=geos, vmin= t_min, vmax= t_max )
        fig.colorbar(mappable=img, label="Brightness Temperature (K)")

        # plot crl path and dots
        # flight path
        os.chdir( crl_path)
        crl_data = xr.open_dataset( crl_name)
        os.chdir( gif_path)
        lon = crl_data.Lon
        lat = crl_data.Lat
        track = sgeom.LineString(zip(lon, lat))
        ax.add_geometries([track], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='g', linewidth=1)
        # legend and title
        sam = mpatches.Rectangle((0, 0), 1, 1, facecolor="g")
        ax.legend([sam], ['CRL Flight Path'])
        plt.title('GOES-16 Clean IR Channel', fontweight='bold', fontsize=15, loc='left')
        plt.title('TC Sam, {}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                  loc='right')
        # ax.tick_params(axis='both',labelsize=200,direction='out',right=False,top=False)
        ax.gridlines(draw_labels=True) # , linewidth=0)

        # dots
        mid_int = ( int( midpoint[11:13]) + int( midpoint[14:16]) / 60 + int(midpoint[17:19] ) / 3600 )
        if mid_int < 10.0:
            mid_int = mid_int + 24.0
        # look at only the decimal places for time ( minutes + seconds)
        for i in range( len( crl_data.time)):
            if ( crl_data.time[i] % 1 <= mid_int % 1 + .01) and ( crl_data.time[i] % 1 >= mid_int % 1 -.01) and (np.rint( crl_data.time[i] ) == np.rint( mid_int)):
                # make star
                long = crl_data.Lon[i].values
                latit = crl_data.Lat[i].values
                ax.scatter( long, latit, s=200, c='g', marker='*', transform=ccrs.PlateCarree() ) # marker = 's'
                break

        # save figure
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/" + str( scan_start_first_dataset.strftime('%m%d%Y')) + "_ir")
        plt.savefig( "goes-image-" + str( goes_ind) + ".png") # str( file_names[ goes_ind][0: -3] ) + '.png', bbox_inches=0)
        os.chdir( gif_path)
        plt.clf()
        print( "Image " + str( goes_ind + 1) + " complete" )



# def goes_visible_gif( flight_path):





def goes_wv_upper_gif( file_names, gif_path, crl_name, crl_path):

    """
    goes_wv_upper_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs

    :param
    :param flight_line: Should this function plot the P-3 flight path above the
                        IR image? Inputs can be either true or false.
    :return: None. The function implicitly saves images to an automatically
             generated folder.
    """

    # make sure that there's a folder in place to hold output images!
    os.chdir( gif_path)
    # use the date from the first image to create a new folder
    first_dataset = xr.open_dataset( file_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    # check if folder already exists. If not, create folder
    os.chdir( "/Users/etmu9498/research/figures/goes-gifs/")
    if not os.path.isdir( str( scan_start_first_dataset.strftime('%m%d%Y')) + "_wv_upper"):
        os.makedirs( "/Users/etmu9498/research/figures/goes-gifs/" + str( scan_start_first_dataset.strftime('%m%d%Y')) + "_wv_upper")
        print( 'New folder created')
    else:
        print( 'Existing folder accessed')
    os.chdir(gif_path)

    # figure out the maximum and minimum brightness temperatures in all datasets
    # to make the colorbar constant between runs
    # just some initial valus that will be automatically rewritten
    t_max = 0
    t_min = 1000
    for goes_ind in range( len( file_names)):
        C = xr.open_dataset( file_names[ goes_ind])
        if np.max( C['CMI_C08'].data) > t_max:
            t_max = np.max( C['CMI_C08'].data )
        if np.min( C['CMI_C08'].data) < t_min:
            t_min = np.min( C['CMI_C08'].data )

    # make images for gif
    for goes_ind in range( len( file_names)):
        # load data for this specific goes_ind
        C = xr.open_dataset( file_names[ goes_ind])

        # Scan's start time, converted to datetime object
        # Scan's end time, converted to datetime object
        # File creation time, convert to datetime object
        # The 't' variable is the scan's midpoint time
        scan_start = datetime.strptime(C.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
        scan_end = datetime.strptime(C.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
        file_created = datetime.strptime(C.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')
        midpoint = str(C['t'].data)[:-8]
        scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

        # access wv data!
        wv_upper = C['CMI_C08'].data
        stacked_wv_upper = np.dstack([ wv_upper, wv_upper, wv_upper])

        print("Starting to make image " + str( goes_ind + 1))

        # make figure
        fig = plt.figure(figsize=(15, 12))
        # these steps use goes ch 2 as a proxy to load things like lat and lon positions
        dat = C.metpy.parse_cf('CMI_C02')
        x = dat.x
        y = dat.y
        geos = dat.metpy.cartopy_crs
        ax = fig.add_subplot(1, 1, 1, projection=geos)

        """
        max_lon = C.geospatial_lat_lon_extent.geospatial_westbound_longitude
        min_lon = C.geospatial_lat_lon_extent.geospatial_eastbound_longitude
        max_lat = C.geospatial_lat_lon_extent.geospatial_northbound_latitude
        min_lat = C.geospatial_lat_lon_extent.geospatial_southbound_latitude
        """

        # plot goes data
        img = ax.imshow( wv_upper, origin='upper',
                  cmap= 'cividis', # plt.cm.get_cmap( 'cividis').reversed(), # "CMRmap",
                  extent=(x.min(), x.max(), y.min(), y.max()),
                  transform=geos, vmin=t_min, vmax=t_max )
        fig.colorbar(mappable=img, label="Brightness Temperature (K)")

        # plot crl path and dots
        # flight path
        os.chdir( crl_path)
        crl_data = xr.open_dataset( crl_name)
        os.chdir( gif_path)
        lon = crl_data.Lon
        lat = crl_data.Lat
        track = sgeom.LineString(zip(lon, lat))
        ax.add_geometries([track], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='g', linewidth=1)
        # legend and title
        sam = mpatches.Rectangle((0, 0), 1, 1, facecolor="g")
        ax.legend([sam], ['CRL Flight Path'])
        plt.title('GOES-16 Upper WV Channel', fontweight='bold', fontsize=15, loc='left')
        plt.title('TC Sam, {}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                  loc='right')
        # ax.tick_params(axis='both',labelsize=200,direction='out',right=False,top=False)
        ax.gridlines(draw_labels=True) # , linewidth=0)

        # dots
        mid_int = ( int( midpoint[11:13]) + int( midpoint[14:16]) / 60 + int(midpoint[17:19] ) / 3600 )
        if mid_int < 10.0:
            mid_int = mid_int + 24.0
        # look at only the decimal places for time ( minutes + seconds)
        for i in range( len( crl_data.time)):
            if ( crl_data.time[i] % 1 <= mid_int % 1 + .01) and ( crl_data.time[i] % 1 >= mid_int % 1 -.01) and (np.rint( crl_data.time[i] ) == np.rint( mid_int)):
                # make star
                long = crl_data.Lon[i].values
                latit = crl_data.Lat[i].values
                ax.scatter( long, latit, s=200, c='g', marker='*', transform=ccrs.PlateCarree() ) # marker = 's'
                break

        # save figure
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/" + str( scan_start_first_dataset.strftime('%m%d%Y')) + "_wv_upper")
        plt.savefig( "goes-image-" + str( goes_ind) + ".png") # str( file_names[ goes_ind][0: -3] ) + '.png', bbox_inches=0)
        os.chdir( gif_path)
        plt.clf()
        print( "Image " + str( goes_ind + 1) + " complete" )





def goes_gif_helper( plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_path, crl_name, crl_path):

    """
    goes_gif_helper is called by the functions above and does most of the work
    for creating the images used in these satellite gifs. It takes specific
    input parameters to customize the output depending on the function calling it
    """

    # make sure that there's a folder in place to hold output images!
    os.chdir( goes_path)
    # use the date from the first image to create a new folder
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    # check if folder already exists. If not, create folder
    os.chdir( "/Users/etmu9498/research/figures/goes-gifs/")
    if not os.path.isdir( output_folder):
        os.makedirs( output_folder)
        print( 'New folder created')
    else:
        print( 'Existing folder accessed')
    os.chdir(goes_path)

    # figure out the maximum and minimum brightness temperatures in all datasets
    # to make the colorbar constant between runs
    # just some initial valus that will be automatically rewritten
    t_max = 0
    t_min = 1000
    for goes_ind in range( len( goes_names)):
        C = xr.open_dataset( goes_names[ goes_ind])
        if np.max( C[ channel_name].data) > t_max:
            t_max = np.max( C[ channel_name].data )
        if np.min( C[ channel_name].data) < t_min:
            t_min = np.min( C[ channel_name].data )

    # make images for gif
    for goes_ind in range( len( goes_names)):
        # load data for this specific goes_ind
        C = xr.open_dataset( goes_names[ goes_ind])

        # Scan's start time, converted to datetime object
        # Scan's end time, converted to datetime object
        # File creation time, convert to datetime object
        # The 't' variable is the scan's midpoint time
        scan_start = datetime.strptime(C.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
        scan_end = datetime.strptime(C.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
        file_created = datetime.strptime(C.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')
        midpoint = str(C['t'].data)[:-8]
        scan_mid = datetime.strptime(midpoint, '%Y-%m-%dT%H:%M:%S.%f')

        # access wv data!
        wv_upper = C[ channel_name].data
        print("Starting to make image " + str( goes_ind + 1))

        # make figure
        fig = plt.figure(figsize=(15, 12))
        # these steps use goes ch 2 as a proxy to load things like lat and lon positions
        dat = C.metpy.parse_cf('CMI_C02')
        x = dat.x
        y = dat.y
        geos = dat.metpy.cartopy_crs
        ax = fig.add_subplot(1, 1, 1, projection=geos)

        """
        failed attempt at trying to set extent of plot
        max_lon = C.geospatial_lat_lon_extent.geospatial_westbound_longitude
        min_lon = C.geospatial_lat_lon_extent.geospatial_eastbound_longitude
        max_lat = C.geospatial_lat_lon_extent.geospatial_northbound_latitude
        min_lat = C.geospatial_lat_lon_extent.geospatial_southbound_latitude
        """

        # plot goes data
        img = ax.imshow( wv_upper, origin='upper',
                  cmap= plot_color, # plt.cm.get_cmap( 'cividis').reversed(), # "CMRmap",
                  extent=(x.min(), x.max(), y.min(), y.max()),
                  transform=geos, vmin=t_min, vmax=t_max )
        fig.colorbar(mappable=img, label="Brightness Temperature (K)")

        # plot crl path and dots
        # flight path
        os.chdir( crl_path)
        crl_data = xr.open_dataset( crl_name)
        os.chdir( goes_path)

        lon = crl_data.Lon
        lat = crl_data.Lat
        track = sgeom.LineString(zip(lon, lat))
        ax.add_geometries([track], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='g', linewidth=1)
        # legend and title
        shape = mpatches.Rectangle((0, 0), 1, 1, facecolor="g")
        ax.legend([ shape], ['CRL Flight Path'])
        plt.title( channel_full_name, fontweight='bold', fontsize=15, loc='left')
        plt.title('TC Sam, {}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                  loc='right')
        # ax.tick_params(axis='both',labelsize=200,direction='out',right=False,top=False)
        ax.gridlines(draw_labels=True) # , linewidth=0)

        # dots
        mid_int = ( int( midpoint[11:13]) + int( midpoint[14:16]) / 60 + int(midpoint[17:19] ) / 3600 )
        if mid_int < 10.0:
            mid_int = mid_int + 24.0
        # look at only the decimal places for time ( minutes + seconds)
        for i in range( len( crl_data.time)):
            if ( crl_data.time[i] % 1 <= mid_int % 1 + .01) and ( crl_data.time[i] % 1 >= mid_int % 1 -.01) and (np.rint( crl_data.time[i] ) == np.rint( mid_int)):
                # make star
                long = crl_data.Lon[i].values
                latit = crl_data.Lat[i].values
                ax.scatter( long, latit, s=200, c='g', marker='*', transform=ccrs.PlateCarree() ) # marker = 's'
                break

        # save figure
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/" + output_folder)
        plt.savefig( "goes-image-" + str( goes_ind) + ".png") # str( file_names[ goes_ind][0: -3] ) + '.png', bbox_inches=0)
        os.chdir( goes_path)
        plt.clf()
        print( "Image " + str( goes_ind + 1) + " complete" )
