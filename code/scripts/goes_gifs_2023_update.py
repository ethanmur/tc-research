"""
Created 1/17/23
Last edited 1/18/23
Description: While the original goes_gif.py code still functions decently, the
file was getting very confusing and clunky, with many locally defined variables
without any explanation. Instead of stressing while trying to simplify and debug
working code, I decided to make a separate, new file that will be vastly simplified!
Especially since I'm moving into analyzing datasets from different years, which
seems like an issue using previous methods.

Edits so far: I really hacked down the goes_helper() and single plot wrapper functions
to make the code much clearer!
"""

# import packages used by these scripts
from datetime import datetime
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy  # noqa: F401
import numpy as np
import xarray as xr
import shapely.geometry as sgeom
import matplotlib.patches as mpatches
from subprocess import run

os.chdir("/Users/etmu9498/research/code/scripts")
import helper_fns
import tc_metadata
import goes_gifs
import make_plots
os.chdir("/Users/etmu9498/research/code/scripts/in-situ-scripts")
import load_in_situ_data

# 1/18 test: this function works!
def load_goes( gif_path, print_files=True):
    return helper_fns.display_data_files( gif_path, 'goes', print_files)

# 1/18 test: this function works!
def print_dates( year):
    if year == 2022:
        goes_folders = [ '0901', '0903', '0904', '0906', '0908', '0916', '0917', '0918', '0920', '0925', '0926', '0927', '1007', '1008']
    elif year == 2021:
        goes_folders = [ '0812am', '0812pm', '0813', '0816', '0817', '0818', '0819', '0820', '0821', '0827', '0926', '0927', '0929']
    else:
        print("Please select either 2022 or 2021 as a valid year for this case!")
        return
    for i in range( len( goes_folders)):
        print( "Dataset "+ str( i) + ": " + str( goes_folders[ i]))


# 1/18 test: this function should work now!
# added functionality includes:
def simple_wrapper( dataset, year=2021, channel='CMI_C13'):
    i = dataset

    if year == 2022:
        goes_folders = [ '0901', '0903', '0904', '0906', '0908', '0916', '0917', '0918', '0920', '0925', '0926', '0927', '1007', '1008']
        
        # crl_path = "/Users/etmu9498/research/data/CRL_data/2022"
        # crl_numbers = [ 0, 2, 4, 5, 6, 7, 8, 47, 51, 58, 60, 10]
        crl_path = "/Users/etmu9498/research/data/cro-all-data-processed/2022"
        crl_numbers = [ 2, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17]

        name_list = [ 'earl', 'earl', 'earl', 'earl', 'earl', 'fiona', 'fiona', 'fiona', 'fiona', 'ian', 'ian', 'ian', 'julia', 'julia']
        extent_list = [ [ -57, -48, 12, 20], [ -65, -57, 15, 23], [ -68, -62, 17, 22], [ -68, -62, 21, 27], [ -68, -62, 25, 31],
                [-62, -48, 8, 20], [-70, -60, 10, 24], [-70, -60, 10, 24], [-75, -64, 14, 27],
                [-85, -75, 12, 20 ], [-85, -75, 12, 24 ], [-89, -79, 18, 26 ], 
                [-90, -70, 0, 20 ], [-90, -70, 0, 20 ]]

        goes_path = "/Users/etmu9498/research/data/goes-satellite/2022/" + goes_folders[ i]

    elif year == 2021:
        goes_folders = [ '0812am', '0812pm', '0813', '0816', '0817', '0818', '0819', '0820', '0821', '0827', '0926', '0927', '0929']

        # load crl data
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        crl_numbers = [ 1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 16, 17, 18]

        in_situ_path = "/Users/etmu9498/research/data/in-situ"
        in_situ_list = make_plots.load_flight_level(in_situ_path, print_files=False)
        in_situ_numbers = [ 4, 5, 6, 7, 8, 9, 12, 14, 16, 17, 31, 32, 34]

        name_list = [ 'fred', 'fred', 'grace', 'grace', 'grace', 'grace', 'henri', 'henri', 'ida', 'sam', 'sam', 'sam']
        extent_list = [ [ -80, -70, 16, 24], [ -80, -70, 16, 24], [-90, -74, 18, 25], # fred zoomed in: [ -77, -73, 20, 24], [ -80, -74, 20, 25],
                [ -76, -68, 15, 20], [ -80, -72, 15, 21], [ -89, -80, 17, 23], [ -95, -85, 18, 25],
                [ -77, -71, 28, 35], [ -75, -67, 34, 40],
                [ -87, -81, 19.5, 25.5],
                [ -55, -46, 12, 18], [ -57, -49, 13, 19], [ -61, -54, 16, 24] ] # the original extent
                # [ -58.5, -57, 19.5, 21] ] # the new, zoom in extent

        # load goes data
        goes_path = "/Users/etmu9498/research/data/goes-satellite/" + goes_folders[ i]


    # outside the loop: set things up
    crl_list = make_plots.load_crl( crl_path, print_files=False)
    goes_names = goes_gifs.load_goes( goes_path, print_files=False)
    crl_name = crl_list[ crl_numbers[ i]]
    extent = extent_list[ i]
    folder = goes_folders[ i]

    # call the helper functions!
    goes_simple( goes_names, goes_path, crl_name, crl_path, folder, year, extent, channel)


# similar to the original goes plotting code below, but it avoids using in situ data and many of the complications found in goes_gif_helper()
def goes_simple( goes_names, goes_path, crl_name, crl_path, folder, year=2021, extent=None, channel_name='CMI_C13'):

    # pick a goes channel for plotting
    if channel_name == 'CMI_C13':
        channel_full_name = "GOES-16 Clean IR Channel"
        plot_color = 'RdYlBu'
    elif channel_name == 'CMI_C08':
        channel_full_name = "GOES-16 Upper Water Vapor Channel"
        plot_color = 'cividis'
    elif channel_name == "CMI_C09":
        channel_full_name = "GOES-16 Mid Water Vapor Channel"
        plot_color = 'plasma'
    elif channel_name == "CMI_C10":
        channel_full_name = "GOES-16 Lower Water Vapor Channel"
        plot_color = 'viridis'
    else:
        print("Please pick a valid GOES channel!")

    # load the first dataset to create the output folder used to save this data!
    os.chdir( goes_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

    if folder == '0812am':
        output_folder = '08122021am_in_situ'
        show_in_situ = False
    elif folder == '0812pm':
        output_folder = '08122021pm_in_situ'
    else:
        output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_in_situ"

    goes_simple_helper( output_folder, channel_name, channel_full_name, goes_names,
        goes_path, crl_name, crl_path, year, extent, plot_color)





# the base function called to create goes plots! This version is a lot simpler than the previous code;
# I got rid of a lot of the confusing options lol.
def goes_simple_helper( output_folder, channel_name, channel_full_name, goes_names,
    goes_path, crl_name, crl_path, year=2021, extent=None, plot_color='RdYlBu', flight_line_color = 'g'):


    # make sure that there's a folder in place to hold output images!
    # check if folder already exists. If not, create folder
    if year == 2022:
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/2022")
    elif year == 2021:
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/")

    if not os.path.isdir( output_folder):
        os.makedirs( output_folder)
        print( 'New folder created: ' + output_folder)
    else:
        print( 'Existing folder accessed: ' + output_folder)
    os.chdir(goes_path)


    # figure out the maximum and minimum brightness temperatures in all datasets
    # to make the colorbar constant between runs!

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
        # print( 'max T: ' + str( t_max))
        # print( 'min T: ' + str( t_min))


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

        # access wv or ir data!
        goes_plot_data = C[ channel_name].data
        print("\nStarting to make image " + str( goes_ind + 1))

        # make figure
        fig = plt.figure(figsize=(15, 12))
        # these steps use goes ch 2 as a proxy to load things like lat and lon positions
        dat = C.metpy.parse_cf('CMI_C02')
        x = dat.x
        y = dat.y
        geos = dat.metpy.cartopy_crs
        ax = fig.add_subplot(1, 1, 1, projection=geos)

        # plot goes data
        img = ax.imshow( goes_plot_data, origin='upper',
                  cmap= plot_color, # plt.cm.get_cmap( 'cividis').reversed(), # "CMRmap",
                  extent=(x.min(), x.max(), y.min(), y.max()),
                  transform=geos, vmin=t_min, vmax=t_max)
        fig.colorbar(mappable=img, label="Brightness Temperature (K)", fraction=0.04, pad=0.05)


        if extent:
            ax.set_extent( extent, crs=ccrs.PlateCarree())


        print( "GOES Image " + str( goes_ind + 1) + " Created")

        # load crl data
        os.chdir( crl_path)
        crl_data = xr.open_dataset( crl_name)
        os.chdir( goes_path)


        # scan time of current goes image
        mid_int = ( int( midpoint[11:13]) + int( midpoint[14:16]) / 60 + int(midpoint[17:19] ) / 3600 )
        print( 'goes time:' + str( mid_int))


        # use a solid line to denote the flight path
        lon = crl_data.Lon.values
        lat = crl_data.Lat.values
        track = sgeom.LineString(zip(lon, lat))
        ax.add_geometries([track], ccrs.PlateCarree(),
                        facecolor='none', edgecolor= flight_line_color, linewidth=1)


        # add dots representing the P-3's location
        lon = crl_data.Lon
        lat = crl_data.Lat
        if mid_int < 10.0:
            mid_int = mid_int + 24.0
        # look at only the decimal places for time ( minutes + seconds)

        star = False
        for i in range( len( crl_data.time)):
            if ( crl_data.time[i] % 1 <= mid_int % 1 + .01) and ( crl_data.time[i] % 1 >= mid_int % 1 -.01) and (np.rint( crl_data.time[i] ) == np.rint( mid_int)):
                # make star
                long = crl_data.Lon[i].values
                latit = crl_data.Lat[i].values
                p3_location = ax.scatter( long, latit, s=125, c= 'k', marker='*', transform=ccrs.PlateCarree() ) # marker = 's'
                star = True
                break

        if star:
            print( "P-3 Location Added")
        else:
            print( 'P-3 location is not within the time frame of this GOES scan')
            p3_location = None


        # add land boundaries to plot?
        # ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
        ax.coastlines()


        # legend and title
        shape = mpatches.Rectangle((0, 0), 1, 1, facecolor=flight_line_color)
        ax.legend([ p3_location, shape], ['P-3 Location', 'CRL Flight Path'])

        plt.title( channel_full_name, fontweight='bold', fontsize=15, loc='left')
        plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                  loc='right')
        # ax.tick_params(axis='both',labelsize=200,direction='out',right=False,top=False)
        ax.gridlines(draw_labels=True) # , linewidth=0)


        # save figure
        if year == 2022:
            os.chdir( "/Users/etmu9498/research/figures/goes-gifs/2022/" + output_folder)
        else:
            os.chdir( "/Users/etmu9498/research/figures/goes-gifs/" + output_folder)

        plt.savefig( "goes-image-" + str( goes_ind) + ".png", dpi=200) # is dpi=200 necessary? # str( file_names[ goes_ind][0: -3] ) + '.png', bbox_inches=0)
        os.chdir( goes_path)
        plt.clf()
        print( "Image " + str( goes_ind + 1) + " complete" )
