# import...
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr
import sys
from scipy.signal import find_peaks
import pandas as pd

os.chdir("/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops
import find_cloud_clusters

# automatically create cloud cluster plots for all datasets given as input!
# much of this code has been taken from "code/scripts/statistics/cloud_clusters.py"
def plot_all_cloud_clusters(tc='all', savefig=False, cluster_type = 'large', cluster_threshold=250, surface_threshold=50): 
    metadata = eyewall_metadata.all_metadata()
    lw = 2
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"

    # use a helper fn to get the relevant years and files
    yearlist, filelist = helper_fns_winter2023.get_crl_datasets( tc=tc)

    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))
    print( 'data saved to: /Users/etmu9498/research/figures/CRL-all-data-processed/"year"-clusters-"cluster type')

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            fs = 20

            # get lidar data for this specific dat
            crl_path = crl_data_root + yearval
            os.chdir( crl_path)
            crl_data = xr.open_dataset( fileval)

            ######
            ## new code: grab the limits for this case
            ######
            date = fileval[7:11]

            # check if this date exists... if not, give it some empty eyewall limits!
            # also check for fred am and pm cases!!
            if date == '0811':
                continue
            if date == '0812':
                if fileval[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                elif fileval[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
            else:
                eyewall_limits = [ ()]

            print("Date: " + date)
            # make plots for proposed eyewall limits!
            for pairi, pairval in enumerate( eyewall_limits):
                print( "Eyewall limits: " + str( pairval))
                # only do this if there are eyewall limits to look at
                if len( pairval) > 1:
                    lefteyewall = pairval[ 0]
                    righteyewall = pairval[ 1]
                # otherwise, give the eyewalls some arbitrary values and don't trim the eyewalls
                else:
                    continue

                ######################
                # plot power ch1
                # also find and plot cloud top heights and clusters!!
                if yearval == '2021':
                    min = -30
                elif yearval == '2022':
                    min = -40
                max=-10

                # preparation
                li = np.argmin( np.abs(crl_data.time.values - lefteyewall))
                ri = np.argmin( np.abs(crl_data.time.values - righteyewall))
                H = crl_data.height
                power = crl_data.P_ch1[ li:ri, :]
                axis = crl_data.time[ li:ri]
                p3_height = crl_data.p3_height[ li:ri]

                # get cloud heights and clusters
                cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = min)
                if cluster_type == 'large':
                    H_clusters, xaxis_clusters = find_cloud_clusters.large_clusters( cloudheights, axis, surface_threshold)
                    H_clusters = H_clusters 
                elif cluster_type == '250m':
                    H_clusters, xaxis_clusters = find_cloud_clusters.clusters_within_250m( cloudheights, axis, cluster_threshold)
                    H_clusters = H_clusters 
                else:
                    print("please enter a valid cloud cluster algorithm name!")
                    return

                ################
                # make the figure!
                plt.figure(figsize=(10, 5))
                helper_fns_winter2023.change_font_sizes(fs, fs)
                plt.title(str(date[0:2] + '/' + date[2:4] + ", Pass " + str(pairi) + " Eye Profile"))


                p = plt.pcolormesh( axis, H, power.transpose(), vmin = min, vmax = max)
                plt.plot( cloudtime, cloudheights, c='r', linewidth=lw, label="cloud tops")
                cbar = plt.colorbar(label="Power Ch. 1 (dBz)", mappable= p, ticks=[-10, -15, -20, -25, -30])
                # Set the number of ticks
                # cbar.locator = ticker.MaxNLocator(nbins=5)
                # cbar.update_ticks()

                plt.legend(loc='upper right')
                plt.ylabel("Height (m)", fontsize=fs)
                plt.xlabel( "Time (Hours, UTC)", fontsize=fs)
                ax=plt.gca()
                ax.set_facecolor('k')
                plt.yticks([1000, 2000, 3000, 4000], fontsize=fs)
                plt.xticks(fontsize=fs)


                # plot mean cloud cluster height lines
                mean_heights = [ None] * len( H_clusters)
                # do this for every cloud cluster
                for i in range( len( H_clusters)):
                    # find the mean value
                    cluster_mean = np.mean( H_clusters[ i])
                    # append the mean value to the list created above
                    # make sure the number of mean heights matches the number of xaxis values for plotting!
                    mean_heights[ i] = [ cluster_mean] * len( H_clusters[ i] )

                for i in range( len( H_clusters)):
                    if i == 0:
                        plt.plot( xaxis_clusters[ i], mean_heights[ i], c='w', linewidth=3, label='Cloud Clusters')
                    else:
                        plt.plot( xaxis_clusters[ i], mean_heights[ i], c='w', linewidth=3)


                if savefig:
                    os.chdir("/Users/etmu9498/research/figures/CRL-all-data-processed/" +yearval + "-clusters-" + cluster_type)
                    plt.savefig(date + '-pass-' + str(pairi) + '.png')