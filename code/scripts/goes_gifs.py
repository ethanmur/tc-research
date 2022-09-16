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

def load_goes( gif_path, print_files=True):
    return helper_fns.display_data_files( gif_path, 'goes', print_files)

def print_dates():
    goes_folders = [ '0812am', '0812pm', '0813', '0816', '0817', '0818', '0819', '0820', '0821', '0827', '0926', '0927', '0929']
    for i in range( len( goes_folders)):
        print( "Dataset "+ str( i) + ": " + str( goes_folders[ i]))

def one_in_situ_plot( dataset):

    # a list of helper information
    goes_folders = [ '0812am', '0812pm', '0813', '0816', '0817', '0818', '0819', '0820', '0821', '0827', '0926', '0927', '0929']

    # load crl data
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    crl_list = make_plots.load_crl( crl_path, print_files=False)
    crl_numbers = [ 1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 16, 17, 18]

    in_situ_path = "/Users/etmu9498/research/data/in-situ"
    in_situ_list = make_plots.load_flight_level(in_situ_path, print_files=False)
    in_situ_numbers = [ 4, 5, 6, 7, 8, 9, 12, 14, 16, 17, 31, 32, 34]

    dataset_list = [ 0, 3,        0, 2, 0, 0, 3, 5] # fix this!! need to edit metadata if arrows are needed
    name_list = [ 'fred', 'fred', 'grace', 'grace', 'grace', 'grace', 'henri', 'henri', 'ida', 'sam', 'sam', 'sam']
    extent_list = [ [ -80, -70, 16, 24], [ -80, -70, 16, 24], [-90, -74, 18, 25], # fred zoomed in: [ -77, -73, 20, 24], [ -80, -74, 20, 25],
            [ -76, -68, 15, 20], [ -80, -72, 15, 21], [ -89, -80, 17, 23], [ -95, -85, 18, 25],
            [ -77, -71, 28, 35], [ -75, -67, 34, 40],
            [ -87, -81, 19.5, 25.5],
            [ -55, -46, 12, 18], [ -57, -49, 13, 19], [ -61, -54, 16, 24] ] # the original extent
            # [ -58.5, -57, 19.5, 21] ] # the new, zoom in extent

    # rename the input dataset as i becuase I'm lazy and don't want to change my old code lol
    i = dataset

    # load goes data
    goes_data_path = "/Users/etmu9498/research/data/goes-satellite/" + goes_folders[ i]
    # goes_data_path = "/Users/etmu9498/research/data/goes-satellite/0829"
    goes_names = goes_gifs.load_goes( goes_data_path, print_files=False)

    crl_name = crl_list[ crl_numbers[ i]]
    in_situ_name = in_situ_list[ in_situ_numbers[ i]]

    # dataset for shear selection
    dataset = dataset_list[ i]
    name = name_list[ i]
    extent = extent_list[ i]
    goes_gifs.goes_in_situ( goes_names, goes_data_path, crl_name, crl_path, extent, in_situ_path, in_situ_name, tcname=name, dataset=dataset, goes_folders= goes_folders[ i])


def all_in_situ_plots():

    # a list of helper information
    goes_folders = [ '0818', '0819', '0820', '0821', '0827', '0926', '0927']

    # load crl data
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    crl_list = make_plots.load_crl( crl_path, print_files=False)
    crl_numbers = [ 7, 8, 9, 11, 12, 16, 17]

    in_situ_path = "/Users/etmu9498/research/data/in-situ"
    in_situ_list = make_plots.load_flight_level(in_situ_path, print_files=False)
    in_situ_numbers = [ 9, 12, 14, 16, 17, 31, 32]

    dataset_list = [ 0, 3, 0, 2, 0, 0, 3]
    name_list = [ 'grace', 'grace', 'henri', 'henri', 'ida', 'sam', 'sam']
    extent_list = [ [ -89, -80, 17, 23], [ -95, -85, 18, 25], [ -77, -71, 28, 35],
            [ -75, -67, 34, 40] , [ -87, -81, 19.5, 25.5], [ -55, -46, 12, 18], [ -57, -49, 13, 19] ]

    for i in range( len( goes_folders)):

        # load goes data
        goes_data_path = "/Users/etmu9498/research/data/goes-satellite/" + goes_folders[ i]
        # goes_data_path = "/Users/etmu9498/research/data/goes-satellite/0829"
        goes_names = goes_gifs.load_goes( goes_data_path, print_files=False)

        crl_name = crl_list[ crl_numbers[ i]]
        in_situ_name = in_situ_list[ in_situ_numbers[ i]]

        # dataset for shear selection
        dataset = dataset_list[ i]
        name = name_list[ i]
        extent = extent_list[ i]
        goes_gifs.goes_in_situ( goes_names, goes_data_path, crl_name, crl_path, extent, in_situ_path, in_situ_name, tcname=name, dataset=dataset)


def goes_in_situ( goes_names, goes_data_path, flight_name, flight_path, extent, p, n, tcname=None, dataset = None, goes_folders= None): # flight_line):
    """
    goes_wv_upper_gif generates images from satellite data (located in a folder
    specified by the user), and saves them automatically to a unique folder.
    The output folder can be found under research/figures/goes-gifs. The clean IR
    band on the goes satellite is used for this analysis.
    """
    channel_name = "CMI_C13" # "CMI_C10" #
    channel_full_name = "GOES-16 Clean IR Channel" # "GOES-16 Lower Water Vapor Channel"#
    os.chdir( goes_data_path)
    first_dataset = xr.open_dataset( goes_names[ 0])
    scan_start_first_dataset = datetime.strptime( first_dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')

    # special cases for fred am and pm flights
    show_in_situ = True

    if goes_folders == '0812am':
        output_folder = '08122021am_in_situ'
        show_in_situ = False
    elif goes_folders == '0812pm':
        output_folder = '08122021pm_in_situ'
    else:
        output_folder = scan_start_first_dataset.strftime('%m%d%Y') + "_in_situ"

    plot_color = 'RdYlBu' # 'Greys'
    flight_line_color = 'g'
    goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name, goes_names,
        goes_data_path, flight_name, flight_path, extent=extent, show_in_situ=show_in_situ, in_situ_data_path=p, in_situ_name=n,
        shear_arrow = False, tcname=tcname, dataset = dataset)



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




def goes_gif_helper( flight_line_color, plot_color, output_folder, channel_name, channel_full_name,
    goes_names, goes_path, crl_name, crl_path, extent=None, show_in_situ=False, in_situ_data_path=None,
    in_situ_name=None, shear_arrow=False, tcname=None, dataset=None):

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


    # load tc data, if applicable
    if tcname and tcname != 'fred':
        metadata = tc_metadata.all_data( tc= tcname)
    else:
        metadata = None



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
        # print( 'max T: ' + str( t_max))
        # print( 'min T: ' + str( t_min))

    # case to make 9/29 data have more contrast
    if output_folder == '09292021_in_situ':
        print( 'special case')
        t_max += - 10
        t_min += 10

        # an extra case accounting for annoying henri errors
        if t_max > 300:
            t_max = 270


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
            scatter_size = 5
        else:
            scatter_size = 2.5

        # scan time of current goes image
        mid_int = ( int( midpoint[11:13]) + int( midpoint[14:16]) / 60 + int(midpoint[17:19] ) / 3600 )
        print( 'goes time:' + str( mid_int))

        # optional code to add shear arrow to plot:

        # also add a dot representing the SHIPS derived tc center, and a p-3 center dot
        add_center = False
        if shear_arrow:

            # load shear data
            shear_mag = metadata[ 'shear_mag'][ dataset]
            shear_dir = metadata[ 'shear_dir'][ dataset] # 90 degrees is from the west

            # print('shear mag: ' + str( shear_mag))
            # print( 'shear direction: ' + str( shear_dir))

            # make the arrow!
            arrow_mag = .2
            # change from the weird shear dir coordinates to sin + cos coordinates
            shear_dir = 90 - shear_dir

            x_tail = 0.1
            y_tail = 0.7
            x_head = x_tail + np.cos( np.deg2rad( shear_dir)) * arrow_mag
            y_head = y_tail + np.sin( np.deg2rad( shear_dir)) * arrow_mag
            arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                 mutation_scale=100,
                                 transform=ax.transAxes, color='k')
            ax.add_patch(arrow)

            # print( 'arrow added')

            shear_label = 'Shear (Black Arrow): ' + str( shear_mag / 10) + ' kt'
            # add text showing shear mag!
            # bbox = {'fc': '0.8', 'pad': 0}
            # props = {'ha': 'center', 'va': 'center', 'bbox': bbox}
            # ax.text( x_tail, y_tail, shear_label, rotation= shear_dir - 10, c='w', transform=ax.transAxes)
            ax.text( .1, .9, shear_label, c='k', transform=ax.transAxes)
            ax.set_ylabel('center / center')
            # (x_tail + x_head) / 2, (y_tail + y_head) / 2

            # print( 'arrow text added')

            # ships_lat = metadata[ 'tc_center_lat'][dataset]
            # ships_lon = metadata[ 'tc_center_lon'][dataset]
            # add_center = True

        # in situ data has been provided case
        if show_in_situ:
            # load data
            in_situ = load_in_situ_data.load_in_situ( in_situ_data_path, in_situ_name, sample_step_size= 1)

            print( 'in situ data loaded')

            # plot flight path line
            lon = in_situ.LONref.values
            lon = [ float( value) for value in lon]
            lat = in_situ.LATref.values
            lat = [ float( value) for value in lat]

            # print( 'lat lons created')

            ws = in_situ["WS.d"].values
            ws = [ float( line) for line in ws]

            # print( 'wind speed created')

            proj = ccrs.PlateCarree()
            img2 = ax.scatter( lon, lat,c = ws, transform= proj, s= scatter_size, marker='o', vmin=0, vmax=60, cmap= 'Greens' )

            # print( 'scatter points added')

            fig.colorbar(mappable=img2, label="Total Wind Speed (m/s)", fraction=0.04, pad=0.05)
            print("Flight Track Added")

            # print( 'colorbar added')

        # otherwise, just use a solid line to denote the flight path
        else:
            lon = crl_data.Lon.values
            lat = crl_data.Lat.values
            track = sgeom.LineString(zip(lon, lat))
            ax.add_geometries([track], ccrs.PlateCarree(),
                            facecolor='none', edgecolor= flight_line_color, linewidth=1)


        # add the ships derived center last for optimal overlap
        # if add_center:
        #     ships_center = ax.scatter( - ships_lon, ships_lat, s=100, c= 'w', marker='p', transform=ccrs.PlateCarree() ) # marker = 's'

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
        # print( 'land boundaries added')


        # legend and title
        shape = mpatches.Rectangle((0, 0), 1, 1, facecolor=flight_line_color)

        if p3_location:
            ax.legend([ p3_location, shape], ['P-3 Location', 'CRL Flight Path'])
        else:
            ax.legend([ shape], [ 'CRL Flight Path'])

        plt.title( channel_full_name, fontweight='bold', fontsize=15, loc='left')
        plt.title('{}'.format(scan_start.strftime('%H:%M UTC %d %B %Y')),
                  loc='right')
        # ax.tick_params(axis='both',labelsize=200,direction='out',right=False,top=False)
        ax.gridlines(draw_labels=True) # , linewidth=0)

        # print( 'saving figure: ')

        # save figure
        os.chdir( "/Users/etmu9498/research/figures/goes-gifs/" + output_folder)
        plt.savefig( "goes-image-" + str( goes_ind) + ".png", dpi=200) # is dpi=200 necessary? # str( file_names[ goes_ind][0: -3] ) + '.png', bbox_inches=0)
        os.chdir( goes_path)
        plt.clf()
        print( "Image " + str( goes_ind + 1) + " complete" )



# this code is meant to show only the flight leg of interest, but it's not working
# so well :(
'''
# only plot in situ data along the specific flight path in question!
# load data
in_situ = load_in_situ_data.load_in_situ( in_situ_data_path, in_situ_name)
in_situ_time = in_situ.time
in_situ_time = [ float( value) for value in in_situ_time]

date = '09-26'
# make a new list of indices for flights only from this day
crl_ind_list = []
for ind_pair in range( len( metadata['crl_range'])):
    if metadata['dates'][ ind_pair] == date:
        crl_ind_list += metadata['crl_range'][ ind_pair]
print( crl_ind_list)
# make a list of times that correspond with this crl index list
crl_time_list = [ crl_data.time[ cur_ind].values for cur_ind in crl_ind_list]
print( crl_time_list)

# find when the goes time falls in between two crl times... only plot
# this in situ data!
lim1 = 0
lim2 = 0
for cur_ind in range( len( crl_time_list) - 1):
    # falling between two crl times
    if mid_int > crl_time_list[ cur_ind] and mid_int <= crl_time_list[ cur_ind + 1]:
        # in situ value closest to 1st time ind
        lim1 = ( np.abs( np.subtract( np.array( in_situ_time), crl_time_list[ cur_ind]))).argmin()
        lim2 = ( np.abs( np.subtract( np.array( in_situ_time), crl_time_list[ cur_ind + 1]))).argmin()
        print( 'limit activated')
    # larger than all crl limits
    elif mid_int > crl_time_list[ cur_ind]:
        lim1 = ( np.abs( np.subtract( np.array( in_situ_time), crl_time_list[ cur_ind]))).argmin()
        lim2 = in_situ_time[ -1]
        print( 'second limit activated')

# print( in_situ_time[ lim1])
# print( in_situ_time[ lim2])
# find the in situ time closest to the goes time
# index = ( np.abs( np.subtract( in_situ_time, mid_int))).argmin()
# print( 'in situ time: ' + str( in_situ_time[ index]))


# plot flight path line
lon = in_situ.LONref[ lim1:lim2]
lon = [ float( value) for value in lon]
lat = in_situ.LATref[ lim1:lim2]
lat = [ float( value) for value in lat]

ws = in_situ["WS.d"].values[ lim1: lim2]
ws = [ float( line) for line in ws]

proj = ccrs.PlateCarree()
img2 = ax.scatter( lon, lat,c = ws, transform= proj, s= scatter_size, marker='o', vmin=0, vmax=60, cmap= 'Greens' )
fig.colorbar(mappable=img2, label="Total Wind Speed (m/s)", fraction=0.04, pad=0.05)

'''
