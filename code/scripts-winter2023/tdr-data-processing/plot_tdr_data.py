import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
os.chdir("/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023

def plot(tc='all'):

	# get a list of datasets to plot
    path = "/Users/etmu9498/research/data/tdr-one-pass/"
    yearlist, namelist, filelist = helper_fns_winter2023.get_tdr_datasets( tc, tdr_data_root=path)

    for yeari, yearval  in enumerate( yearlist):
        for tci, tcval in enumerate( namelist[yeari]):
            for filei, fileval in enumerate( filelist[yeari][tci][0]):
                # print("year " + yearval + ", tc " + tcval + ", file " + fileval)
                
                # get data
                os.chdir(path + "/" + yearval + "/" + tcval)
                tdr_data = xr.open_dataset(fileval)

                # make plot
                plt.figure(figsize=(10, 3))
                color_map = plt.cm.get_cmap( "RdYlBu").reversed()
                xaxis = tdr_data.distance
                sz = 16
                helper_fns_winter2023.change_font_sizes(sz, sz)
                # plot data
                # get rid of nans and resize array to get rid of overlapping data
                # also, no need to use .transpose() because that was already done when making the datasets!
                reflectivity = tdr_data.REFLECTIVITY[ :, 0:len( xaxis)]
                # refl = tdr_data.REFLECTIVITY[:, lat_no_nan_ind] # another way?
                plt.pcolormesh( xaxis, tdr_data.height, reflectivity, cmap = color_map, vmin = -10, vmax = 50 )
                cbar = plt.colorbar( )
                cbar.ax.set_ylabel( ylabel="TDR Reflectivity (dBZ)", fontsize=sz)
                cbar.ax.tick_params(labelsize=sz)
                plt.ylabel( 'Height (Km)')
                plt.xlabel( 'Radial Distance (Km)')
                plt.title('TDR Data: TC ' + tcval + ", " + fileval[-16:-3])

                # save the plot
                # see if there's already a tc name folder availible
                os.chdir("/Users/etmu9498/research/figures/tdr-one-pass/")
                output_folder = tcval
                if not os.path.isdir( output_folder):
                    os.makedirs( output_folder)
                    print( 'New folder created: ' + output_folder)

                savedir = "/Users/etmu9498/research/figures/tdr-one-pass/" + tcval
                os.chdir( savedir)
                date = fileval[-16:-3]
                filename = 'tdr_combined_' +tcval+ '_' +date+ '.png'
                plt.savefig( filename, dpi=200, bbox_inches='tight')
