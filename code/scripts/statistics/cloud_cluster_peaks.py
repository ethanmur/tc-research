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
def clusters_all_peaks( tc='all', no_eyewalls=True, cluster_threshold=.250, print_stats=True):
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
            clusters_one_tc_eye( new_crl_path, new_crl_name, eyewall_dists, title, cluster_threshold)

            save_dir = "/Users/etmu9498/research/figures/prob-dist-results/clusters-peaks"

            os.chdir( save_dir)
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )

            # os.chdir( save_path)
            # plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )

    return




# use the new crl datasets; working with distance coordinates will be much easier!
def clusters_one_tc_eye( new_crl_path, new_crl_name, eyewall_dists, title, cluster_threshold):

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

    Hpeaks, xpeaks = cloud_cluster_algorithms.find_peaks_clusters( H, xaxis, cluster_threshold)

    #######################
    ## Statistics
    #######################



    ############ plotting ##################
    # plot the crl backscattered power figure for helpful testing
    # plt.figure( figsize=(18, 5))
    fig, a0 = plt.subplots(1, 1, figsize=(14, 5), facecolor='w')
    helper_fns.change_font_sizes(small=14, medium=16 )
    plt.sca( a0)
    make_plots_new_heights.plot_power_ch1( new_crl_path, new_crl_name, 'in-situ-dist')

    # plot cloud top heights from algorithm
    a0.scatter( xaxis, H, c= 'r', s=11, marker='s')
    # a0.plot( xaxis, H, c=  'r', linewidth=2, label= 'Cloud Top Height')

    # plot mean cloud cluster height lines
    a0.scatter( xpeaks, Hpeaks, c='w', s=60, marker='*', label='Cluster Peaks')

    a0.set_xlabel( 'Distance from TC Center ( Km)')
    a0.set_xlim( eyewall_dists)
    a0.legend(loc='upper right', framealpha=1, facecolor='.85', bbox_to_anchor=(1.01, 1.15))
    a0.set_title( title, loc='left')
    a0.set_ylim( [ 0, 3.8])
