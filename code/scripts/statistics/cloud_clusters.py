import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns # for making prettier pdfs
import warnings
import pandas as pd

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import cloud_height
import tc_metadata
import helper_fns

os.chdir(  "/Users/etmu9498/research/code/scripts/statistics")
import cloud_cluster_algorithms

# wrapper function used to run and save plots for all tc eye passes
def clusters_all_tc_eyes( tc='all', no_eyewalls=True, cluster_threshold=.250, print_stats=True, surface_cases=True, new_algorithm=True):
    warnings.filterwarnings("ignore")

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # look at a specific tc, or every tc!
    for tcname in tcname_list:
        # load data
        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return

        # print some helpful notices to the user
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')

        # an empty list that will hold all the height datasets from within the eye
        # all_heights = []

        # make a probability distribution function for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):
            # load data from new sources
            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)
            # new_tdr_path = metadata[ 'new_tdr_path']
            new_crl_path = metadata[ 'um_crl_path']

            # new code for eyewall vs no eyewall in situ distances!
            if no_eyewalls:
                eyewall_dists = metadata[ 'eyewall_dists_no_eyewalls'][ dataset]
                # save_path = "/Users/etmu9498/research/figures/prob-dist-results/clusters-no-eyewalls/"
            else:
                # save_path = "/Users/etmu9498/research/figures/prob-dist-results/clusters-in-situ/"
                eyewall_dists = metadata[ 'in_situ_eyewall_dists'][ dataset]

            title = ( "CRL Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # print out the current dataset for the user!
            print( "TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset])

            # make and save the figure for this pass!
            clusters_one_tc_eye( new_crl_path, new_crl_name, eyewall_dists, title, cluster_threshold, print_stats, surface_cases, new_algorithm)

            if new_algorithm:
                save_dir = "/Users/etmu9498/research/figures/prob-dist-results/clusters-larger"
            else:
                save_dir = "/Users/etmu9498/research/figures/prob-dist-results/clusters-no-eyewalls"

            os.chdir( save_dir)
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )

            # os.chdir( save_path)
            # plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )

    return

'''
            # code to look at heights (or eventually clusters) of clouds for all tc eye passes for one tc
            # implement this later


            # save height values from this run
            all_heights = all_heights + np.ndarray.tolist( H_dataset)

        # return statistics for all cases averaged together!
        all_heights = np.array( all_heights)
        print( "Number of data points:  " + str( len( all_heights)))
        print( "Height value range:     " + str( all_heights.min()) + " km to " + str( all_heights.max()) + " km")
        print( "Height value mean:      " + str( np.mean( all_heights)) + " km")
        print( "Height value median:    " + str( np.median( all_heights)) + " km\n")

        # make a histogram representing the total cloud layer for this tc!
        plt.figure( figsize=(8, 6), facecolor='w')
        sns.distplot( all_heights, bins=25, hist=True, vertical=True, color='y')

        plt.xlabel( 'Probability of a Given Height Occurence')
        plt.ylabel( 'Height from Surface (Km)')
        plt.title( "Histogram of All Cloud Heights for TC " + metadata['tc_name'])
        plt.ylim( [0, 3.8])
        plt.xlim( [0, 2])
        plt.grid('on')

        # save the histogram
        os.chdir( save_path)
        plt.savefig( metadata['tc_name'].casefold() + "-total.png", bbox_inches='tight', dpi=300 )

'''

# use the new crl datasets; working with distance coordinates will be much easier!
def clusters_one_tc_eye( new_crl_path, new_crl_name, eyewall_dists, title, cluster_threshold, print_stats, surface_cases, new_algorithm):

    # load data
    os.chdir( new_crl_path)
    crl_data = xr.open_dataset( new_crl_name)
    xaxis_data = crl_data.in_situ_distance.values

    # find the indices and values in the crl distance dataset closest to the eyewall_dists limits
    i1, x1 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 0])
    i2, x2 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 1])

    # find cloud top heights for values within the specified eye distance range
    H, xaxis = cloud_height.find_cloud_heights(new_crl_name, -30, i1, i2, xaxis ='in-situ-dist', crl_path=new_crl_path, new_heights=True)

    #################################
    ## new code for finding cloud clusters: use helper algorithm
    #################################

    if new_algorithm:
        print("Using new algorithm")
        H_clusters, xaxis_clusters = cloud_cluster_algorithms.large_clusters( H, xaxis, cluster_threshold)
    else:
        print("Using new algorithm")
        H_clusters, xaxis_clusters = cloud_cluster_algorithms.clusters_within_250m( H, xaxis, cluster_threshold)

    #######################
    ## Statistics
    #######################

    #######################
    ## cluster mean heights
    #######################

    # find mean height values for each cloud cluster: save them in the same length format as the clusters for plotting!
    # initialize a mean height array
    # the array is initialized with empty values and filled in over time
    mean_heights = [ None] * len( H_clusters)

    # a smaller array holding only one value per cloud cluster.
    # meant for statistical calculations, not plotting
    mean_height_vals = [ None] * len( H_clusters)

    # do this for every cloud cluster
    for i in range( len( H_clusters)):
        # find the mean value
        cluster_mean = np.mean( H_clusters[ i])
        # append the mean value to the list created above
        # make sure the number of mean heights matches the number of xaxis values for plotting!
        mean_heights[ i] = [ cluster_mean] * len( H_clusters[ i] )

        mean_height_vals[ i] = cluster_mean


    #######################
    ## cluster mean widths
    #######################

    # find the widths of each cluster
    mean_width_vals = [ None] * len( xaxis_clusters)
    # do this for every cloud cluster
    for i in range( len( xaxis_clusters)):

        # save the current and next clusters in temp variables
        cluster = xaxis_clusters[ i]
        # nextcluster = xaxis_clusters[ i+1]

        # cases for finding the delta distance
        # two positive distances: simply take the difference
        if cluster[ 0] >= 0 and cluster[ -1] >= 0:
            mean_width_vals[ i] = abs( cluster[ -1] - cluster[ 0])

        # one positive distance: flip the negative sign with a -
        elif cluster[ 0] >= 0 and cluster[ -1] < 0:
            mean_width_vals[ i] = abs( cluster[ 0] - cluster[ -1])

        # one positive distance: flip the negative sign with a -
        elif cluster[ 0] < 0 and cluster[ -1] >= 0:
            mean_width_vals[ i] = abs( cluster[ -1] - cluster[ 0])

        # two negative distances: flip both the negative signs with a -
        elif cluster[ 0] < 0 and cluster[ -1] < 0:
            mean_width_vals[ i] = abs( cluster[ -1] - cluster[ 0])

        else:
            print( 'error when finding mean cell widths')


    ############ plotting ##################
    # plot the crl backscattered power figure for helpful testing
    # plt.figure( figsize=(18, 5))
    fig, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 1]}, figsize=(12, 9), facecolor='w')
    helper_fns.change_font_sizes(small=14, medium=16 )
    plt.sca( a0)
    make_plots_new_heights.plot_power_ch1( new_crl_path, new_crl_name, 'in-situ-dist')

    # plot cloud top heights from algorithm
    a0.scatter( xaxis, H, c= 'r', s=11, marker='s')
    a0.plot( xaxis, H, c=  'r', linewidth=2, label= 'Cloud Top Height')



    # special case accounting for surface clouds
    if surface_cases:

        #######################
        ## separate sea spray + surface clusters from clouds!
        #######################

        xaxis_surf_clusters = np.array( xaxis_clusters) [ np.where( np.array( mean_height_vals) < 0.05)[0] ]
        xaxis_cloud_clusters = np.array( xaxis_clusters) [ np.where( np.array( mean_height_vals) >= 0.05)[0] ]
        h_surf_clusters = np.array( H_clusters) [ np.where( np.array( mean_height_vals) < 0.05)[0] ]
        h_cloud_clusters = np.array( H_clusters) [ np.where( np.array( mean_height_vals) >= 0.05)[0] ]
        mean_surf_h = np.array( mean_heights) [ np.where( np.array( mean_height_vals) < 0.05)[0] ]
        mean_cloud_h = np.array( mean_heights) [ np.where( np.array( mean_height_vals) >= 0.05)[0] ]
        mean_surf_width = np.array( mean_width_vals) [ np.where( np.array( mean_height_vals) < 0.05)[0] ]
        mean_cloud_width = np.array( mean_width_vals) [ np.where( np.array( mean_height_vals) >= 0.05)[0] ]

        '''
        # testing
        lists = [ xaxis_surf_clusters, xaxis_cloud_clusters, h_surf_clusters, h_cloud_clusters, mean_surf, mean_cloud]
        for i in lists:
            print( len( i))
        '''

        # plot mean cloud cluster height lines: distinguish between surface and normal clouds!
        for i in range( len( h_surf_clusters)):
            plt.plot( xaxis_surf_clusters[ i], np.array( mean_surf_h[ i]) + .05, c='w', linewidth=3)

        # plot mean cloud cluster height lines: distinguish between surface and normal clouds!
        for i in range( len( h_cloud_clusters)):
            plt.plot( xaxis_cloud_clusters[ i], np.array( mean_cloud_h[ i]), c='b', linewidth=3)

        # make an empty line to add cloud cluster label to legend
        plt.plot( -2, -2, c='w', linewidth=3, label='Surface Cluster Mean Height')
        plt.plot( -2, -2, c='b', linewidth=3, label='Cloud Cluster Mean Height')


    # normal case: assume surface layer can be a cloud cluster
    else:

        # plot mean cloud cluster height lines
        for i in range( len( H_clusters)):
            plt.plot( xaxis_clusters[ i], mean_heights[ i], c='w', linewidth=3)
        # make an empty line to add cloud cluster label to legend
        plt.plot( -2, -2, c='w', linewidth=3, label='Cloud Cluster Mean Height')


    a0.set_xlabel( 'Distance from TC Center ( Km)')
    a0.set_xlim( eyewall_dists)
    a0.legend(loc='upper right', framealpha=1, facecolor='.85', bbox_to_anchor=(1.01, 1.15))
    a0.set_title( title, loc='left')
    a0.set_ylim( [ 0, 3.8])


    # view simple histogram of cloud cluster widths and heights!
    nbins = 25 # number of bins

    plt.sca( a1)


    if surface_cases:

        # maybe use this code earlier to speed things up!!
        # i did :) the code to efficiently find surface vs cloud heights is above now

        # remove clusters of only 1 point from analysis
        surf_widths = mean_surf_width [ np.where( mean_surf_width > 0.0)[0] ]
        cloud_widths = mean_cloud_width [ np.where( mean_cloud_width > 0.0)[0] ]

        surf_widths = surf_widths.tolist()
        cloud_widths = cloud_widths.tolist()

        # print( surf_widths)
        # print( cloud_widths)
        # print( type( surf_widths))


        # turn width arrays into a pandas dataframe for plotting!
        # not working ://///
        '''
        dataset = {
            'name' : ['Surface Widths', 'Cloud Widths'],
            'Surface Widths' : surf_widths,
            'Cloud Widths' : cloud_widths
            }

        # creating a Dataframe object
        df = pd.DataFrame( dataset)
        '''

        # sns.set(style='white')
        # df.set_index('data').plot(kind='bar', stacked=True, color=['steelblue', 'red'])

        sns.histplot( x= cloud_widths, element='bars', edgecolor='k', linewidth=3, color='b', kde=True, bins=nbins, stat="count") # color = 'g', edgecolor='k', linewidth=2,

    # include zero cases in statistics
    else:
        # remove clusters of only 1 point from analysis
        nonzero_widths = np.array( mean_width_vals) [ np.where( np.array( mean_width_vals) > 0.0)[0] ]
        sns.histplot( x= nonzero_widths, kde=True, bins=nbins, stat="count", edgecolor='k', linewidth=2, color = 'g')

    a1.set_ylabel( 'Cell Count at each Cluster Width', fontsize=16)
    a1.set_xlabel( 'Cluster Width (Km)', fontsize=16)
    # a1.set_xlim( [-0.2, 30.0])
    # a1.set_xlim( [0, 2])
    a1.grid(True)
    a1.tick_params(axis='both', which='major', labelsize=16)

    if print_stats:
        print( "heights: values and count")
        print( len( mean_height_vals))
        print( mean_height_vals)
        print( "widths: values and count")
        print( len( mean_width_vals))
        print( mean_width_vals)

        print( "widths without 0 values: count")
        print( len( np.where( np.array( mean_width_vals) > 0.0)[0]))
        print( '\n\n')


    '''
    # old code for plotting surface vs tall clouds!

    xaxis_surf_clusters = []
    xaxis_cloud_clusters = []
    h_surf_clusters = []
    h_cloud_clusters = []
    mean_surf = []
    mean_cloud = []

    # a clear surface is considered when the "clouds" are 50m tall or less- this could
    # just be p-3 flight errors or aerosols or tall waves?
    for i in range( len( mean_height_vals)):

        # surface cluster case
        if mean_height_vals[ i] < .05:
            xaxis_surf_clusters.append( xaxis_clusters[ i])
            h_surf_clusters.append( H_clusters[ i])
            mean_surf.append( mean_heights[ i])

        # above 50 m case: cloud clusters
        else:
            xaxis_cloud_clusters.append( xaxis_clusters[ i])
            h_cloud_clusters.append( H_clusters[ i])
            mean_cloud.append( mean_heights[ i])
    '''
