import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import seaborn as sns # for making prettier pdfs
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import cloud_height
import tc_metadata
import helper_fns


def pdf_all_tc_eyes( tc='all', no_eyewalls=False):
    warnings.filterwarnings("ignore")

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['fred', 'grace', 'henri', 'ida', 'sam']
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
        all_heights = []

        # make a probability distribution function for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):
            # load data from new sources
            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)
            # new_tdr_path = metadata[ 'new_tdr_path']
            new_crl_path = metadata[ 'um_crl_path']

            # new code for eyewall vs no eyewall in situ distances!
            if no_eyewalls:
                eyewall_dists = metadata[ 'eyewall_dists_no_eyewalls'][ dataset]
                save_path = "/Users/etmu9498/research/figures/prob-dist-results/pdfs-v1-no-eyewalls/"
            else:
                save_path = "/Users/etmu9498/research/figures/prob-dist-results/pdfs-v1-in-situ/"
                eyewall_dists = metadata[ 'in_situ_eyewall_dists'][ dataset]

            shear_quad = metadata['shear_quads'][dataset]

            title = ( "CRL Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # print out the current dataset for the user!
            print( "TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset])

            H_dataset = pdf_one_tc_eye( new_crl_path, new_crl_name, eyewall_dists, title, shear_quad)

            os.chdir( save_path)
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )

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

    return

# use the new crl datasets; working with distance coordinates will be much easier!
def pdf_one_tc_eye( new_crl_path, new_crl_name, eyewall_dists, title, shear_quad):

    # load data
    os.chdir( new_crl_path)
    crl_data = xr.open_dataset( new_crl_name)
    xaxis_data = crl_data.in_situ_distance.values

    # find the indices and values in the crl distance dataset closest to the eyewall_dists limits
    i1, x1 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 0])
    i2, x2 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 1])

    # print( "i1 and x1: " + str( i1) + " " + str(x1))
    # print( "i2 and x2: " + str( i2) + " " + str(x2))

    # print( 'i1 and x1: ' + str( i1) + ' ' + str( x1))
    # print( 'i2 and x2: ' + str( i2) + ' ' + str( x2))

    # find cloud top heights for values within the specified eye distance range
    H, xaxis_value = cloud_height.find_cloud_heights(new_crl_name, -30, i1, i2, xaxis ='in-situ-dist', crl_path=new_crl_path, new_heights=True)

    # plot the crl backscattered power figure for helpful testing
    # plt.figure( figsize=(18, 5))
    fig, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [1, 2, .75]}, figsize=(12, 18), facecolor='w')
    helper_fns.change_font_sizes(small=14, medium=16 )

    plt.sca( a0)
    make_plots_new_heights.plot_power_ch1( new_crl_path, new_crl_name, 'in-situ-dist')
    a0.scatter( xaxis_value, H, c= 'r', s=11, marker='s')
    a0.plot( xaxis_value, H, c=  'r', linewidth=2, label= 'Cloud Top Height')

    a0.set_xlabel( 'Distance from TC Center ( Km)')
    a0.set_xlim( eyewall_dists)
    a0.legend(loc='lower left')
    a0.set_title( title)
    a0.set_ylim( [ 0, 3.8])

    # add shear info as text!

    label0 = '0'
    label1 = '1'

    if shear_quad[0] == 'UL':
        label0 = 'Upshear Left'
    elif shear_quad[0] == 'UR':
        label0 = 'Upshear Right'
    elif shear_quad[0] == 'DL':
        label0 = 'Downshear Left'
    elif shear_quad[0] == 'DR':
        label0 = 'Downshear Right'

    if shear_quad[1] == 'UL':
        label1 = 'Upshear Left'
    elif shear_quad[1] == 'UR':
        label1 = 'Upshear Right'
    elif shear_quad[1] == 'DL':
        label1 = 'Downshear Left'
    elif shear_quad[1] == 'DR':
        label1 = 'Downshear Right'

    # only add first label if axis contains negative values ( to the left of 0 km)
    if eyewall_dists[0] < 0: # 3.25
        a0.text( .025, .9, label0, color='k', fontsize=12,
                bbox={'facecolor': 'w', 'alpha': 0.85, 'pad': 10},
                verticalalignment='bottom', horizontalalignment='left',
                transform=a0.transAxes)

    a0.text( .975, .9, label1, color='k', fontsize=12,
            bbox={'facecolor': 'w', 'alpha': 0.85, 'pad': 10},
            verticalalignment='bottom', horizontalalignment='right',
            transform=a0.transAxes)


    # view simple histogram of cloud top heights!
    nbins = 25 # number of bins

    # plt.hist( H, nbins, density=True, facecolor='c', histtype = 'bar') # 'bar', 'barstacked', 'step', 'stepfilled'

    plt.sca( a1)
    sns.distplot( H, bins=nbins, hist=True, vertical=True, color='y', norm_hist=False)

    # print( H)

    # plt.axvline( x=H.mean(), c='g', label="Mean height value", linewidth=3)

    a1.set_xlabel( 'Probability of a Given Height Occurence')
    a1.set_ylabel( 'Height from Surface (Km)')
    # plt.title( "TC Sam, 09/26/22, Eye 1, Cloud Top Height Histogram")
    a1.set_ylim( [-0.2, 3.8])
    a1.set_xlim( [0, 2])
    a1.grid(True)

    # a clear surface is considered when the "clouds" are 50m tall or less- this could
    # just be p-3 flight errors or aerosols or tall waves?
    if len( H) > 0:
        cloud_frac = ( len( np.where( H < .05 )[0]) / len( H) )
        cloud_percent = cloud_frac * 100

        # save some basic stats from the height dataset to image!
        a2.text( 0, .95, ("Number of data points:  " + str( i2 - i1)))
        a2.text( 0, .8, ("Number of clear air data points:  " + str( len( np.where( H < .05 )[0]) )))
        a2.text( 0, .6, ("Height value range:     " + str( H.min()) + " km to " + str( H.max()) + " km"))
        a2.text( 0, .4, ("Height value mean:      " + str( H.mean()) + " km"))
        a2.text( 0, .2, ("Height value median:    " + str( np.median( H)) + " km\n"))
        a2.text( 0, .1, ("Clear air percentage: " + str( cloud_percent) + " %"))
        a2.set_axis_off()
    else:
        a2.text( 0, .95, "Error when calculating cloud heights :(")


    return H






def number_of_layers( tc='all'):
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

            # new code for in situ eyewall distances!
            eyewall_dists = metadata[ 'in_situ_eyewall_dists'][ dataset]

            title = ( "CRL Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # print out the current dataset for the user!
            print( "TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset])

            H_dataset = number_of_layers_one_eye( new_crl_path, new_crl_name, eyewall_dists, title)

            os.chdir( "/Users/etmu9498/research/figures/pdfs-number-of-layers/")
            plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=300 )

            '''
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
        os.chdir( "/Users/etmu9498/research/figures/pdfs-v1-in-situ/")
        plt.savefig( metadata['tc_name'].casefold() + "-total.png", bbox_inches='tight', dpi=300 )
        '''
    return

# use the new crl datasets; working with distance coordinates will be much easier!
def number_of_layers_one_eye( new_crl_path, new_crl_name, eyewall_dists, title):

    # load data
    os.chdir( new_crl_path)
    crl_data = xr.open_dataset( new_crl_name)
    xaxis_data = crl_data.in_situ_distance.values

    # find the indices and values in the crl distance dataset closest to the eyewall_dists limits
    i1, x1 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 0])
    i2, x2 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 1])

    # find cloud top heights and layer counts for values within the specified eye distance range
    H, xaxis_value, cloud_counts = cloud_height.find_multi_cloud_heights(new_crl_name, -30, i1, i2, xaxis ='in-situ-dist', crl_path=new_crl_path)

    layers = range( len( cloud_counts))
    print( "Number of cloud layers: " + str( layers))
    print( "Count: " + str( cloud_counts))

    # plot the crl backscattered power figure for helpful testing
    fig, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 2]}, figsize=(12, 12), facecolor='w')
    helper_fns.change_font_sizes(small=14, medium=16 )
    plt.sca( a0)
    make_plots_new_heights.plot_power_ch1( new_crl_path, new_crl_name, 'in-situ-dist')

    # trying to make hollow markers but this won't work for some reason??!
    a0.scatter( xaxis_value, H, c= 'r', s=15) # fillstyle='none') # facecolors='none', edgecolors='r') # marker='s')
    # a0.plot( xaxis_value, H, markerfacecolor='none', ms=15, markeredgecolor='red')

    a0.set_xlabel( 'Distance from TC Center ( Km)')
    a0.set_xlim( eyewall_dists)
    a0.legend(loc='lower left')
    a0.set_title( title)
    a0.set_ylim( [ 0, 3.2])

    # view simple histogram of cloud top layer counts!
    # nbins = 10 # number of bins
    plt.sca( a1)
    # sns.histplot( data=np.array(cloud_counts), x=np.array(layers), binwidth=1, color='c', stat='count') # bins=nbins
    # sns.barplot( x=layers, y=cloud_counts, color='c') # binwidth=1, , stat='count') # bins=nbins
    # plt.axvline( x=H.mean(), c='g', label="Mean height value", linewidth=3)
    plt.bar( x=layers, height=cloud_counts, color='c')

    a1.set_xlabel( 'Number of cloud layers')
    a1.set_ylabel( 'Count')
    # a1.set_yscale( 'log')
    a1.grid(True)

    return H





# this function finds statistics on cloud heights depending on the tc intensity!
# it outputs figures summarizing these values for each intensity category
def cloud_height_vs_intensity( csu_poster_case=False, no_eyewalls=False):
    warnings.filterwarnings("ignore")

    # empty lists that will hold all the height datasets for each intensity category
    td_heights, td_cases = [], 0
    ts_heights, ts_cases = [], 0
    wh_heights, wh_cases = [], 0
    sh_heights, sh_cases = [], 0
    # cycle through every dataset
    tcname_list = [ 'fred', 'grace', 'henri', 'ida', 'sam']
    for tcname in tcname_list:
        # load data
        metadata = tc_metadata.all_data( tcname)
        if metadata == 'selected TC name is not yet implemented':
            print( metadata)
            return
        # print some helpful notices to the user
        print( "\nTC " + metadata['tc_name'])
        print( 'Number of crl files: ' + str( len( metadata['xlims'] ))+ '\n')

        # get cloud top height information for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):
            # print out the current dataset for the user
            print( "TC " + metadata['tc_name'] + ", " + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset])

            # setup
            # load data from new sources

            # use clipped data for 4 ind cases
            # normal case
            # if metadata['crl_range'][dataset] == 2:
            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)
            new_crl_path = metadata[ 'um_crl_path']
            # elif metadata['crl_range'][dataset] == 4:
            #     new_crl_path = metadata[ 'um_crl_path']


            # new code for eyewall vs no eyewall in situ distances!
            if no_eyewalls:
                eyewall_dists = metadata[ 'eyewall_dists_no_eyewalls'][ dataset]
            else:
                eyewall_dists = metadata[ 'in_situ_eyewall_dists'][ dataset]

            title = ( "CRL Data, TC " + metadata['tc_name'] + ", "
                    + metadata['dates'][ dataset] + ", Eye Pass " + metadata['eye_pass'][ dataset] )

            # This bit of code is from the pdf_one_tc_eye() function
            os.chdir( new_crl_path)
            crl_data = xr.open_dataset( new_crl_name)
            xaxis_data = crl_data.in_situ_distance.values
            i1, x1 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 0])
            i2, x2 = helper_fns.closest_val( xaxis_data, eyewall_dists[ 1])

            # find cloud top heights for values within the specified eye distance range
            H, xaxis_value = cloud_height.find_cloud_heights(new_crl_name, -30, i1, i2, xaxis ='in-situ-dist', crl_path=new_crl_path, new_heights=True)


            # save height values from this run
            cat = metadata['intensity_cat'][dataset] # tc intensity category
            # figure out where to put the height data depending on the tc intensity!
            if cat == 'td':
                td_heights += np.ndarray.tolist( H)
                # print( 'td')
                td_cases += 1
            elif cat == 'ts':
                # print( 'ts')
                ts_heights += np.ndarray.tolist( H)
                ts_cases += 1
            elif cat == 'wh':
                # print( 'wh')
                wh_heights += np.ndarray.tolist( H)
                wh_cases += 1
            elif cat == 'sh':
                # print( 'sh')
                sh_heights += np.ndarray.tolist( H)
                sh_cases += 1

    # create histograms for every TC intensity
    # put things in lists for easier looping
    fig_title = [ 'tropical-depressions', 'tropical-storms', 'weak-hurricanes', 'strong-hurricanes']
    # nicely formatted titles for plots, as opposed to saving
    fig_title_nice = [ 'Tropical Depressions', 'Tropical Storms', 'Weak Hurricanes', 'Strong Hurricanes']

    cases = [ td_cases, ts_cases, wh_cases, sh_cases]
    heights = [ td_heights, ts_heights, wh_heights, sh_heights]
    # loop through each case
    for i in range( 4):
        height = np.array( heights[ i])

        # make pretty figures for csu poster conference!
        if csu_poster_case:

            # remove 0 km heights from figure!
            plot_height = height[ np.where( height > .05)[0] ]

            # no cases for this category (only applies to td's because Fred case hasn't been added yet)
            if len( height) != 0:

                helper_fns.change_font_sizes(small=24, medium=24 )
                fig, a0 = plt.subplots(1, 1, figsize=( 8.5, 7) )

                # set the figure background (not plot background!) to transparent!!
                fig.patch.set_facecolor('blue')
                fig.patch.set_alpha(0)

                # create a histogram
                plt.sca( a0)

                # sns.displot( x=height, kind='kde')
                sns.set_theme(style="white", palette=None)

                # the line below was probs causing the annoying lack of plot space!!!
                # sns.set_context(rc = {'patch.linewidth': 0.0})
                sns.histplot( y= plot_height, kde=True, binwidth=.175, edgecolor='k', linewidth=2, color = 'g') # color = 'g', element='poly') # hist_kws=dict(edgecolor="black", linewidth=2))
                a0.set_xlabel( 'Count for Each Cloud Height')
                a0.set_ylabel( 'Height from Surface (Km)')
                a0.set_title( "Cloud Height Distribution for " + fig_title_nice[ i], fontsize=24)
                a0.set_ylim( [-.2, 3.75])
                a0.set_xlim( [0, 200])
                a0.grid(True)

                # save the histogram
                os.chdir( "/Users/etmu9498/research/figures/csu-poster/")
                if no_eyewalls:
                    plt.savefig( fig_title[i] + "-no-eyewalls.png", bbox_inches='tight', dpi=500, transparent=False )
                else:
                    plt.savefig( fig_title[i] + "-new.png", bbox_inches='tight', dpi=500, transparent=False )
                print( fig_title[i] + ' figure created.')

            else:
                print( "The tropical depression case hasn't yet been added!")
                continue

        # normal figure creation case
        else:
            # no cases for this category (only applies to td's because Fred case hasn't been added yet)
            if len( height) != 0:
                fig, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, figsize=(10, 12), facecolor='w')

                helper_fns.change_font_sizes(small=18, medium=18 )

                # create a histogram
                plt.sca( a0)

                # sns.displot( x=height, kind='kde')
                sns.set_theme(style="white", palette=None)
                sns.set_context(rc = {'patch.linewidth': 0.0})
                sns.histplot( y=height, kde=True, bins=25, color = 'y')
                a0.set_xlabel( 'Count of Given Cloud Height')
                a0.set_ylabel( 'Height from Surface (Km)')
                a0.set_title( "Histogram of All Cloud Heights for " + fig_title_nice[ i])
                a0.set_ylim( [-.2, 3.5])
                a0.set_xlim( [0, 250])
                a0.grid(True)
                a0.axhline( y=height.mean(), c='g', label="Mean height value", linewidth=3)
                a0.axhline( y=np.median(height), c='b', label="Median height value", linewidth=3)
                a0.legend(loc='upper right')

                '''
                # a0.set_xlim( [0, 1.3])
                sns.distplot( height, bins=25, hist=True, vertical=True, color='y')
                a0.set_xlabel( 'Probability of a Given Height Occurence')
                '''

                # print out some useful statistics
                plot_height = height[ np.where( height > .05)[0] ]
                clear_points  = len( height[ np.where( height < .05)[0] ])

                a1.text( 0, 1, ("Number of data points:  " + str( len( height))))
                a1.text( 0, .8, ("Number of clear air points: " + str( clear_points )))
                # a1.text( 0, .8, ("Height value range:     " + str( height.min()) + " km to " + str( height.max()) + " km"))
                a1.text( 0, .6, ("Height value mean (no clear air): " + str( plot_height.mean()) + " km"))
                # a1.text( 0, .3, ("Height value median:    " + str( np.median( height)) + " km\n"))
                a1.text( 0, .3, ("clear air fraction:    " + str( 100 * ( clear_points / len(height) ))))
                a1.text( 0, .2, ("Number of eye passes: " + str( cases[ i])))
                a1.set_axis_off()
                # save the histogram
                os.chdir( "/Users/etmu9498/research/figures/prob-dist-results/pdfs-intensity/")
                if no_eyewalls:
                    plt.savefig( fig_title[i] + "-no-eyewalls.png", bbox_inches='tight', dpi=300 )
                else:
                    plt.savefig( fig_title[i] + "-new.png", bbox_inches='tight', dpi=300 )
                print( fig_title[i] + ' figure created.')

            else:
                print( "The tropical depression case hasn't yet been added!")
                continue
