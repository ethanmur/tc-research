# import...
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr
import sys
import scipy.signal as signal
import pandas as pd
import random

os.chdir("/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops
import find_cloud_clusters

# automatically create cloud cluster plots for all datasets given as input!
# much of this code has been taken from "code/scripts/statistics/cloud_clusters.py"
def plot_all_cloud_clusters(tc='all', savefig=False, makefig=False, cluster_type = 'large', cluster_threshold=250, surface_threshold=50, min_size=False, low_pass_filt=False, combine_clusters=False): 
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

    # save flight / pass info and cloud cluster sizes here!
    df_peaks = pd.DataFrame( )
    # save values in temporary lists here: will be added to df_peaks once filled
    df_yearlist = []
    df_datelist = []
    df_passlist = []
    df_clusterlist = []
    df_meanlist = []

    print(filelist)

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            fs = 20

            # get lidar data for this specific date
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

                # optional code to apply a low pass filter to the height data! maybe better results?
                if low_pass_filt:
                    ### Apply a Butterworth filter (recursive filter) to the data
                    N = 9 ## order
                    frequency_cutoff = .20 # 0.15 # remove f's greater than .05 year^-1 (20 years)
                    Wn = frequency_cutoff*2 ## scalar given the critical frequency (all higher frequencies are removed)

                    ## Note: Wn is normalized from 0 to 1, where 1 is the Nyquist frequency, pi radians/sample. 
                    ## Note: (Wn is thus in half-cycles / sample.) 
                    b, a = signal.butter(N, Wn)

                    # for plotting later
                    orig_heights = cloudheights.copy()

                    cloudheights = signal.filtfilt(b, a, cloudheights) ## one filter forward, one filter backward - you are filtering twice
                    # make sure filtered clouds don't dip below 0m
                    cloudheights[np.where(cloudheights < 0.0)[0]] = 0.0


                if cluster_type == 'large':
                    H_clusters, xaxis_clusters = find_cloud_clusters.large_clusters( cloudheights, axis, surface_threshold)
                elif cluster_type == '250m':
                    H_clusters, xaxis_clusters = find_cloud_clusters.clusters_within_250m( cloudheights, axis, cluster_threshold, surface_threshold)
                elif cluster_type == '250m-improved':
                    H_clusters, xaxis_clusters = find_cloud_clusters.clusters_within_threshold_updated( cloudheights, axis, cluster_threshold, surface_threshold, min_size=min_size, combine_clusters=combine_clusters)
                
                    # print(H_clusters)
                    # print(xaxis_clusters)

                elif cluster_type == '250m-multi':
                    cloudheights, cloudtime, cloud_counts = find_cloud_tops.find_multi_cloud_heights( H, power, axis, p3_height, cutoff_power = min)
                    H_clusters = find_cloud_clusters.multi_clusters_within_250m( cloudheights, axis)
                else:
                    print("please enter a valid cloud cluster algorithm name!")
                    return

                # find mean cloud cluster heights
                mean_heights = [ None] * len( H_clusters)
                cluster_means = [None] * len( H_clusters)
                # do this for every cloud cluster
                for i in range( len( H_clusters)):
                    # print("cluster " + str(i))
                    # print(H_clusters[ i])

                    # find the mean value
                    # or just plot clouds all at the top of the figure! with a small nudge
                    # cluster_mean = 3000 + 500 * random.random()
                    cluster_mean = np.nanmean( np.array( H_clusters[ i]))
                    cluster_means[i] = cluster_mean / 1000.
                    
                    # append the mean value to the list created above
                    # make sure the number of mean heights matches the number of xaxis values for plotting!
                    mean_heights[ i] = [ cluster_mean] * len( H_clusters[ i] )
                   
                # find cloud cluster horizontal sizes
                cluster_size = []
                for i in range( len( H_clusters)):
                    cluster_size.append( ((2 * 130) / 1000) * len( H_clusters[i]))

                ################
                # make the figure!
                if makefig:
                    plt.figure(figsize=(10, 10))

                    plt.subplot(211)
                    helper_fns_winter2023.change_font_sizes(fs, fs)
                    plt.title(str(date[0:2] + '/' + date[2:4] + ", Pass " + str(pairi) + " Eye Profile"))


                    p = plt.pcolormesh( axis, H, power.transpose(), vmin = min, vmax = max)

                    if low_pass_filt:
                        plt.plot( cloudtime, orig_heights, c='c', linewidth=lw, label="original clouds")
                    
                    plt.plot( cloudtime, cloudheights, c='r', linewidth=lw, label="cloud tops")

                    cbar = plt.colorbar(label="Power Ch. 1 (dBz)", mappable= p, ticks=[-10, -15, -20, -25, -30])
                    # Set the number of ticks
                    # cbar.locator = ticker.MaxNLocator(nbins=5)
                    # cbar.update_ticks()

                    # plt.legend(loc='upper right')
                    plt.ylabel("Height (m)", fontsize=fs)
                    plt.xlabel( "Time (Hours, UTC)", fontsize=fs)
                    ax=plt.gca()
                    ax.set_facecolor('k')
                    plt.yticks([1000, 2000, 3000, 4000], fontsize=fs)
                    plt.xticks(fontsize=fs)
                    
                    for i in range( len( H_clusters)):
                        # if i == 0:
                        #     plt.scatter( xaxis_clusters[ i], mean_heights[ i], c='w', s=5, label='Cloud Clusters', marker='s')
                        # else:

                        # print(len(xaxis_clusters[i]))
                        # print(len(mean_heights[i]))

                        print('case ' + str(i))
                        print( xaxis_clusters[ i][0] - .05)
                        print( mean_heights[i][0] - 200)
                        print(xaxis_clusters[ i][-1] - xaxis_clusters[ i][0] + .05)
                        plt.Rectangle( (xaxis_clusters[ i][0] - .05, mean_heights[i][0] - 200), xaxis_clusters[ i][-1] - xaxis_clusters[ i][0] + .05, 500, fc='w',ec="w")
                        plt.plot( xaxis_clusters[ i], mean_heights[ i], c='w', linewidth=3) # , marker='s')



                    # add cloud cluster size histogram!
                    plt.subplot(212)

                    # same size as cloud height bins!
                    bw1 = (2*130) / 1000

                    # special case for clusters larger than 20 km
                    if len(cluster_size) <= 0:
                        xinc1 = np.arange(0, 20, bw1)
                    elif np.nanmax( np.array(cluster_size)) > 20:
                        print("large cluster case")
                        xinc1 = np.arange(0, np.nanmax(np.array(cluster_size)) + 2.0, bw1)
                    else:
                        xinc1 = np.arange(0, 20, bw1)
                    
                    hx=np.histogram( cluster_size,xinc1)
                    plt.rcParams["figure.figsize"] = [5,5]
                    plt.bar(hx[1][:-1],hx[0],edgecolor = 'k', color = 'forestgreen', width = bw1, linewidth = 1)
                    plt.ylabel('Count')
                    plt.xlabel('Cloud Cluster Size (Km)')

                # save info for dataframe later!
                df_yearlist.append( yearval)
                df_datelist.append( date)
                df_passlist.append( pairi)
                df_clusterlist.append( cluster_size)
                df_meanlist.append( cluster_means)

                if savefig:
                    os.chdir("/Users/etmu9498/research/figures/CRL-all-data-processed/" +yearval + "-clusters-" + cluster_type)
                    plt.savefig(date + '-pass-' + str(pairi) + '.png')

    # return the cluster dataframe
    df_peaks['year'] = df_yearlist
    df_peaks['date'] = df_datelist
    df_peaks['pass'] = df_passlist
    df_peaks['cluster sizes (km)'] = df_clusterlist
    df_peaks['cluster mean heights (km)'] = df_meanlist

    return df_peaks