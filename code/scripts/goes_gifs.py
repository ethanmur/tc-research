"""
Created 6/1/22
Last edited 6/9/22
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
import plot_in_situ

def load_goes( gif_path, print_files=True):
    return helper_fns.display_data_files( gif_path, 'goes', print_files)


def goes_in_situ( goes_names, goes_data_path, flight_name, flight_path, extent, p, n): # flight_line):
    """
    goes_wv_upper_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The clean IR
    band on the goes satellite is used for this analysis.
    """
    channel_name = "CMI_C09" # "CMI_C13"
    channel_full_name = "GOES-16 Mid Water Vapor Channel"# "GOES-16 Clean IR Channel"
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_in_situ"
    plot_color = 'Greys'
    flight_line_color = 'g'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_data_path, flight_name, flight_path, extent=extent, show_in_situ=True, in_situ_data_path=p, in_situ_name=n)



def goes_ir_gif( goes_names, goes_data_path, crl_name, crl_path): # flight_line):
    """
    goes_wv_upper_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The clean IR
    band on the goes satellite is used for this analysis.
    """
    channel_name = "CMI_C13"
    channel_full_name = "GOES-16 Clean IR Channel"
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_ir"
    plot_color = 'Greys'
    flight_line_color = 'k'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_data_path, crl_name, crl_path)


# def goes_visible_gif( flight_path):

def goes_wv_upper_gif( goes_names, goes_data_path, crl_name, crl_path):
    """
    goes_wv_upper_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The upper
    atmosphere water vapor observing band on the goes satellite is used for this
    analysis.
    """
    channel_name = "CMI_C08"
    channel_full_name = "GOES-16 Upper Water Vapor Channel"
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_wv_upper"
    plot_color = cmocean.cm.haline # 'cividis'
    flight_line_color = 'w'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_data_path, crl_name, crl_path)


def goes_wv_mid_gif( goes_names, goes_data_path, crl_name, crl_path):
    """
    goes_wv_mid_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The middle
    atmosphere water vapor observing band on the goes satellite is used for this
    analysis.
    """
    channel_name = "CMI_C09"
    channel_full_name = "GOES-16 Mid Water Vapor Channel"
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_wv_mid"
    plot_color = 'viridis'
    flight_line_color = 'w'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_data_path, crl_name, crl_path)


def goes_wv_lower_gif( goes_names, goes_data_path, crl_name, crl_path):
    """
    goes_wv_lower_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The lower
    atmosphere water vapor observing band on the goes satellite is used for this
    analysis.
    """
    channel_name = "CMI_C10"
    channel_full_name = "GOES-16 Lower Water Vapor Channel"
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_wv_lower"
    plot_color = 'plasma'
    flight_line_color = 'w'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_data_path, crl_name, crl_path)


def collocation_test( goes_names, goes_data_path, crl_name, crl_path):
    """
    goes_wv_lower_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The lower
    atmosphere water vapor observing band on the goes satellite is used for this
    analysis.
    """
    channel_name = "CMI_C10"
    channel_full_name = "GOES-16 Lower Water Vapor Channel"
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    output_folder = "collocation-test"
    plot_color = 'plasma'
    flight_line_color = 'w'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_data_path, crl_name, crl_path)




def goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names, goes_path, crl_name, crl_path, extent=None, show_in_situ=False, in_situ_data_path=None, in_situ_name=None):

    """
    goes_gif_helper is called by the functions above and does most of the work
    for creating the images used in these satellite gifs. It takes specific
    input parameters to customize the output depending on the function calling it
    """

    # make sure that there's a folder in place to hold output images!

    # check if folder already exists. If not, create folder
    os.chdir( "/Users/etmu9498/research/figures/goes-gifs/")
    if not os.path.isdir( output_folder):
        os.makedirs( output_folder)
        print( 'New folder created: ' + output_folder)
    else:
        print( 'Existing folder accessed: ' + output_folder)
    os.chdir(goes_path)

    # figure out the maximum and minimum brightness temperatures in all datasets
    # to make the colorbar constant between runs

    # just some initial values that will be automatically rewritten. 0 is tiny for a max
    # (literally absolute zero lol) so it must be rewritten by the first dataset.
    # then, the goal is to find the max and min out of all datasets, and use those values
    # as constant colorbar min and maxes
    t_max = 0
    t_min = 10000
    for goes_ind in range( len( goes_names)):
        C = xr.open_dataset( goes_names[ goes_ind])

        if np.nanmax( C[ channel_name].data) > t_max:
            t_max = np.nanmax( C[ channel_name].data )
        if np.nanmin( C[ channel_name].data) < t_min:
            t_min = np.nanmin( C[ channel_name].data )

    print("Number of GOES Images to be Created: " + str( len( goes_names)))

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
        print("\nStarting to make image " + str( goes_ind + 1))

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
                  transform=geos, vmin=t_min, vmax=t_max)
        fig.colorbar(mappable=img, label="Brightness Temperature (K)", fraction=0.04, pad=0.05)

        print( "GOES Image " + str( goes_ind + 1) + " Created")

        # load crl data
        os.chdir( crl_path)
        crl_data = xr.open_dataset( crl_name)
        os.chdir( goes_path)

        # optional function call to crop data
        if extent:
            ax.set_extent( extent, crs=ccrs.PlateCarree())
            # make this size bigger when zooming in
            scatter_size = 15
        else:
            scatter_size = 5

        # in situ data has been provided case
        if show_in_situ:
            # load data
            in_situ = plot_in_situ.load_in_situ( in_situ_data_path, in_situ_name)

            # plot flight path line
            lon = in_situ.LONref
            lon = [ float( value) for value in lon]
            lat = in_situ.LATref
            lat = [ float( value) for value in lat]

            ws = in_situ["WS.d"].values
            ws = [ float( line) for line in ws]

            proj = ccrs.PlateCarree()
            img2 = ax.scatter( lon, lat,c = ws, transform= proj, s= scatter_size, marker='o', vmin=0, vmax=60, cmap= 'Greens' )
            fig.colorbar(mappable=img2, label="Total Wind Speed (m/s)", fraction=0.04, pad=0.05)

        # otherwise, just use a solid line to denote the flight path
        else:
            track = sgeom.LineString(zip(lon, lat))
            ax.add_geometries([track], ccrs.PlateCarree(),
                            facecolor='none', edgecolor= flight_line_color, linewidth=1)

        print("Flight Track Added")

        # add dots representing the P-3's location
        lon = crl_data.Lon
        lat = crl_data.Lat
        mid_int = ( int( midpoint[11:13]) + int( midpoint[14:16]) / 60 + int(midpoint[17:19] ) / 3600 )
        if mid_int < 10.0:
            mid_int = mid_int + 24.0
        # look at only the decimal places for time ( minutes + seconds)
        for i in range( len( crl_data.time)):
            if ( crl_data.time[i] % 1 <= mid_int % 1 + .01) and ( crl_data.time[i] % 1 >= mid_int % 1 -.01) and (np.rint( crl_data.time[i] ) == np.rint( mid_int)):
                # make star
                long = crl_data.Lon[i].values
                latit = crl_data.Lat[i].values
                ax.scatter( long, latit, s=200, c= flight_line_color, marker='*', transform=ccrs.PlateCarree() ) # marker = 's'
                break

        print( "P-3 Location Added")


        # legend and title
        shape = mpatches.Rectangle((0, 0), 1, 1, facecolor=flight_line_color)
        ax.legend([ shape], ['CRL Flight Path'])
        plt.title( channel_full_name, fontweight='bold', fontsize=15, loc='left')
        plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                  loc='right')
        # ax.tick_params(axis='both',labelsize=200,direction='out',right=False,top=False)
        ax.gridlines(draw_labels=True) # , linewidth=0)


        # save figure
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/" + output_folder)
        plt.savefig( "goes-image-" + str( goes_ind) + ".png", dpi=200) # is dpi=200 necessary? # str( file_names[ goes_ind][0: -3] ) + '.png', bbox_inches=0)
        os.chdir( goes_path)
        plt.clf()
        print( "Image " + str( goes_ind + 1) + " complete" )
