import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import matplotlib.ticker as ticker


import xarray as xr
import seaborn as sns # for making prettier pdfs
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import cloud_height
import tc_metadata
import helper_fns








# this function finds statistics on cloud heights depending on the tc intensity!
# it outputs figures summarizing these values for each intensity category
def cloud_height_vs_intensity( no_eyewalls=True, binwidth=1.0, percent_plot=False, density_plot=False, poster_case=False, lw=2, manual_bins=False, smoothwidth=10):
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

    colors = [ 'b', 'k', 'y', 'g']

    cases = [ td_cases, ts_cases, wh_cases, sh_cases]
    heights = [ td_heights, ts_heights, wh_heights, sh_heights]



    # make figure before loop
    fig, a0 = plt.subplots(1, 1, figsize=( 10, 10) )
    smallfont=18
    mediumfont=20
    largefont=24

    # set the figure background (not plot background!) to transparent!!
    fig.patch.set_facecolor('blue')
    fig.patch.set_alpha(0)

    # create a histogram
    plt.sca( a0)



    # print(len( heights))
    # print( len( heights[ 0]))

    # find the total number of height records
    total_heights = 0
    for i in range( 4):
        total_heights += len( heights[ i])


    # loop through each case
    for i in range( 4):
        height = np.array( heights[ i])

        # remove 0 km heights from figure!
        plot_height = height[ np.where( height > .05)[0] ]

        # no cases for this category (only applies to td's because Fred case hasn't been added yet)
        if len( height) != 0:

            # just some semantic definition changes... not important
            width = binwidth

            # special case for making a nice figure for the esss poster conference!
            if poster_case:
                # this step is to have all integrals of the curves to sum to 1!
                # not just to have each curve sum to 1
                # a0.xaxis.set_major_formatter( PercentFormatter( 1 * scale_factor, symbol=None, decimals=2))
                # ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format( plot_height * scale_factor))
                # a0.xaxis.set_major_formatter(ticks)
                # cut=0 tells the function to not plot curves above the highest and lowest measured values
                # try using clip as well! might be better. see docs


                # new method based on zhien's code: create manual bins for stats, don't use seaborn
                # don't use a gaussian curve to interpolate data; just use regular bins
                if manual_bins:
                    savedir = "/Users/etmu9498/research/figures/prob-dist-results/intensity-same-plot-new/"

                    height_bin=np.arange(0, 4.51, width)

                    # print( height_bin)

                    # calculate mean height variation within the bin range
                    mean_height=[]
                    mean_height_count=[]
                    mean_height_prob = []

                    mean_height_prob_trimmed = []
                    height_bin_trimmed = []

                    # a list that will hold max smoothed heights from all loops
                    # used for proper x axis scaling
                    max_heights = []

                    # print( len( plot_height))

                    # do this for every height level determined by the manually inputed bin width
                    for newi in height_bin:

                        # find the points that fall within this height bin for this step
                        res=np.where(np.logical_and( plot_height >= newi - width / 2., plot_height <= newi + width / 2. ))

                        # print( plot_height[ res])

                        mean_height.append( np.mean( plot_height[ res]))
                        mean_height_count.append( len( res[0]))
                        # use height, not plot_height, to include heights below 50m in caclulations
                        # this line accounts for clear air fraction when scaling the curves!
                        mean_height_prob.append( len( res[0] ) / len( height))



                    # print( len( height))
                    # print( len( plot_height))
                    # print( mean_height)
                    # print( mean_height_count)

                    # clip cloud height counts and bins so that only full bins are plotted: ignore empty bins


                    # smooth data before plotting to eliminate noise
                    box_pts = smoothwidth
                    box = np.ones(box_pts)/box_pts
                    prob_smooth = np.convolve( mean_height_prob, box, mode='same')

                    # a0.plot( mean_height_prob, height_bin, color=colors[i], linewidth=lw/2, alpha=.5)
                    a0.plot( prob_smooth * 100, height_bin, color=colors[i], label=fig_title_nice[i], linewidth=lw)


                    a0.set_xlabel('Cloud Height Probability (%)', fontsize=mediumfont)
                    # a0.set_xlabel('Cloud Height Probability (Sums to 1)', fontsize=mediumfont)

                    # use the highest smoothed value to effectively scale the x axis!!
                    max_prob = np.nanmax( prob_smooth)
                    max_heights.append( max_prob)

                    # print( 'max prob: ' + str( max_prob))


                    # make and save separate plots for each intensity! in the same output folder
                    fig1, a1 = plt.subplots(1, 1, figsize=( 10, 10) )
                    # set the figure background (not plot background!) to transparent!!
                    fig1.patch.set_facecolor('blue')
                    fig1.patch.set_alpha(0)
                    plt.sca( a1)
                    a1.plot( np.array( mean_height_prob) * 100, height_bin, color=colors[i], label='Raw Data: Bins = ' + str( width * 1000) + " m", linewidth=lw/2, alpha=.7)
                    a1.plot( prob_smooth * 100, height_bin, color=colors[i], label='Smoothed Data', linewidth=lw)

                    # a1.set_xlabel('Cloud Height Probability (Sums to 1)', fontsize=mediumfont)
                    a1.set_xlabel('Cloud Height Probability (%)', fontsize=mediumfont)

                    # sns.set_theme(style="white", palette=None)
                    a1.set_ylabel( 'Height from Surface (Km)', fontsize=mediumfont)
                    a1.set_ylim( [-.2, 4.5])

                    # a0.set_xlim( [ 0, 1.05])
                    # a1.set_xlim( [ 0, .02])
                    a1.set_xlim( [ 0, max_prob * 2 * 100])


                    a1.grid(True)
                    a1.set_title( "Cloud Height Distributions for " + fig_title_nice[i], fontsize=largefont)
                    a1.tick_params(axis='both', which='major', labelsize=smallfont)
                    a1.legend(fontsize=smallfont, loc="upper right")

                    os.chdir( savedir)
                    plt.savefig( fig_title[i] + '-' + str( int( binwidth*1000)) + "m.png", bbox_inches='tight', dpi=500, transparent=False )



                    plt.sca( a0)


                # original method
                else:
                    savedir = "/Users/etmu9498/research/figures/prob-dist-results/intensity-same-plot/"
                    sns.kdeplot( y= plot_height, bw_adjust= width, linewidth=lw, color = colors[i], label=fig_title_nice[i], cut=0) #  clip=[0, 3.6]) # bw_adjust=1 -> curve found in histplot

                    a0.set_xlabel('Cloud Height Probability Density', fontsize=mediumfont)

                    a0.set_xlim( [ 0, 1.05])



                # sns.set_theme(style="white", palette=None)
                a0.set_ylabel( 'Height from Surface (Km)', fontsize=mediumfont)
                a0.set_ylim( [-.2, 4.5])

                a0.grid(True)
                a0.set_title( "Cloud Height Distributions by TC Intensity", fontsize=largefont)
                a0.tick_params(axis='both', which='major', labelsize=smallfont)
                plt.legend(fontsize=smallfont, loc="upper right")


                # save figs at the end of the function!
                continue



            # show percentages on the xaxis
            if percent_plot:
                a0.xaxis.set_major_formatter( PercentFormatter(1 / width))  # show axis such that 1/binwidth corresponds to 100%
                a0.set_xlabel(f'Cloud Height Percents for Bin Width = {width} km', fontsize=mediumfont)
                sns.kdeplot( y= plot_height, bw_adjust= width, linewidth=lw, color = colors[i], label=fig_title_nice[i]) # bw_adjust=1 -> curve found in histplot

            # show densities on the xaxis (area under curve sums to 1)
            elif density_plot:
                sns.kdeplot( y= plot_height, bw_adjust= width, linewidth=lw, color = colors[i], label=fig_title_nice[i]) # bw_adjust=1 -> curve found in histplot
                a0.set_xlabel(f'Cloud Height Density for Bin Width = {width} km', fontsize=mediumfont)

            # show height probabilities (percent / 100)
            else:
                a0.xaxis.set_major_formatter( PercentFormatter(100 / (width), symbol=None, decimals=2))  # show axis such that 1/binwidth corresponds to 100%
                # plot_height = plot_height / width
                a0.set_xlabel(f'Cloud Height Probability for Bin Width = {width} km', fontsize=mediumfont)
                sns.kdeplot( y=  plot_height, bw_adjust= width, linewidth=lw, color = colors[i], label=fig_title_nice[i]) # bw_adjust=1 -> curve found in histplot


    if poster_case:
        # make one last auto axis adjustment for new, manual plots
        if manual_bins:
            # find the highest max value, and scale the x axis with that
            total_max = np.max( max_heights)
            plt.sca( a0)

            a0.set_xlim( [ 0, total_max * 1.5 * 100])

            # 25 m case
            # a0.set_xlim( [ 0, 1.6])
            # plt.xticks([ 0, .5, 1.0, 1.5])

            # 100 m case
            # a0.set_xlim( [0, 7])


            # print( total_max)

        # save the histogram
        os.chdir( savedir)
        plt.savefig( "one-fig-pdfs-" + str( int( binwidth*1000)) + "m.png", bbox_inches='tight', dpi=500, transparent=False )

        return



    # sns.displot( x=height, kind='kde')
    sns.set_theme(style="white", palette=None)

    # a0.set_xlabel( 'Percentage of CRL Returns at a Given Height')
    a0.set_ylabel( 'Height from Surface (Km)', fontsize=mediumfont)
    # a0.set_ylim( [-.2, 3.75])
    a0.set_ylim( [-.2, 4.5])
    # a0.set_xlim( [0, 200])
    # a0.set_xlim( [-.05, 2])
    a0.grid(True)
    a0.set_title( "Cloud Height Distributions by TC Intensity", fontsize=largefont)
    a0.tick_params(axis='both', which='major', labelsize=smallfont)
    plt.legend(fontsize=smallfont, loc="upper right")

    # saving options
    if percent_plot:
        # save the histogram
        os.chdir( "/Users/etmu9498/research/figures/prob-dist-results/intensity-same-plot/")
        if no_eyewalls:
            plt.savefig( "all-dists-no-eyewalls-" + str( binwidth) + "-percent.png", bbox_inches='tight', dpi=500, transparent=False )
        else:
            plt.savefig( "all-dists-" + str( binwidth) + "-percent.png", bbox_inches='tight', dpi=500, transparent=False )


    elif density_plot:
        # save the histogram
        os.chdir( "/Users/etmu9498/research/figures/prob-dist-results/intensity-same-plot/")
        if no_eyewalls:
            plt.savefig( "all-dists-no-eyewalls-" + str( binwidth) + "-density.png", bbox_inches='tight', dpi=500, transparent=False )
        else:
            plt.savefig( "all-dists-" + str( binwidth) + "-density.png", bbox_inches='tight', dpi=500, transparent=False )

    else:
        # save the histogram
        os.chdir( "/Users/etmu9498/research/figures/prob-dist-results/intensity-same-plot/")
        if no_eyewalls:
            plt.savefig( "all-dists-no-eyewalls-" + str( binwidth) + "-probability.png", bbox_inches='tight', dpi=500, transparent=False )
        else:
            plt.savefig( "all-dists-" + str( binwidth) + "-probability.png", bbox_inches='tight', dpi=500, transparent=False )






# pretty much the same code as above, but with separate distribution profile panels based on intensity
def by_intensity_stacked( binwidth=1.0):
    warnings.filterwarnings("ignore")

    # empty lists that will hold all the height datasets for each intensity category
    td_heights, td_cases = [], 0
    ts_heights, ts_cases = [], 0
    wh_heights, wh_cases = [], 0
    sh_heights, sh_cases = [], 0

    # cycle through every dataset
    tcname_list = ['fred', 'grace', 'henri', 'ida', 'sam']
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

            tdr_name, new_crl_name = tc_metadata.choose_new_data( tcname, dataset)
            new_crl_path = metadata[ 'um_crl_path']

            eyewall_dists = metadata[ 'eyewall_dists_no_eyewalls'][ dataset]

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

    colors = [ 'b', 'k', 'g', 'c']

    cases = [ td_cases, ts_cases, wh_cases, sh_cases]
    heights = [ td_heights, ts_heights, wh_heights, sh_heights]


    # make figure before loop
    fig = plt.figure( figsize=( 6.5, 12.5) )
    smallfont=14
    mediumfont=16
    largefont=20

    # set the figure background (not plot background!) to transparent!!
    fig.patch.set_facecolor('blue')
    fig.patch.set_alpha(0)

    subplot_list = [ 411, 412, 413, 414]

    # loop through each case
    for i in range( 4):

        plt.subplot( subplot_list[ i])

        height = np.array( heights[ i])

        # remove 0 km heights from figure!
        plot_height = height[ np.where( height > .05)[0] ]

        width = binwidth

        # show height probabilities (percent / 100)
        a0 = plt.gca()

        a0.xaxis.set_major_formatter( PercentFormatter(100 / (width), symbol=None, decimals=2))  # show axis such that 1/binwidth corresponds to 100%
        # plot_height = plot_height / width
        sns.kdeplot( y=  plot_height, bw_adjust= width, linewidth=2, color = colors[i], label=fig_title_nice[i]) # bw_adjust=1 -> curve found in histplot

        if i == 0:
            a0.set_title( "Cloud Height Distributions by TC Intensity", fontsize=largefont)
            plt.legend(fontsize=smallfont, loc="lower right")

        else:
            plt.legend(fontsize=smallfont, loc="upper right")

        if i < 3:
            # a0.set_xticklabels([])
            a0.set_xticklabels([])
        else:
            a0.set_xlabel(f'Cloud Height Probability', fontsize=mediumfont)

        a0.set_xlim([0, 1.2])


        sns.set_theme(style="white", palette=None)

        # a0.set_xlabel( 'Percentage of CRL Returns at a Given Height')
        # a0.set_ylabel( 'Height from Surface (Km)', fontsize=mediumfont)
        a0.tick_params(axis='both', which='major', labelsize=16)

        a0.set_ylim( [-.2, 3.75])
        a0.grid(True)

    # a0.tick_params(axis='both', which='major', labelsize=smallfont)

    # manually set y axis label text
    a0.text( -.15, 8, "Height from Surface (Km)", c='k', fontsize=mediumfont, rotation=90,
           verticalalignment='center', horizontalalignment='center')


    plt.subplots_adjust(wspace=0, hspace=0)

    os.chdir( "/Users/etmu9498/research/figures/esss-poster/")
    plt.savefig( "stacked-pdfs.png", bbox_inches='tight', dpi=500, transparent=False )






# this function finds statistics on cloud heights depending on the tc intensity!
# it outputs figures summarizing these values for each intensity category
def cloud_height_vs_intensity_testing( no_eyewalls=True):
    warnings.filterwarnings("ignore")

    # empty lists that will hold all the height datasets for each intensity category
    td_heights, td_cases = [], 0
    ts_heights, ts_cases = [], 0
    wh_heights, wh_cases = [], 0
    sh_heights, sh_cases = [], 0

    # cycle through every dataset
    tcname_list = ['grace', 'henri', 'ida', 'sam']
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

    colors = [ 'b', 'k', 'g', 'c']

    cases = [ td_cases, ts_cases, wh_cases, sh_cases]
    heights = [ td_heights, ts_heights, wh_heights, sh_heights]



    # loop through each case
    for i in range( 4):
        height = np.array( heights[ i])

        # remove 0 km heights from figure!
        plot_height = height[ np.where( height > .05)[0] ]

        # no cases for this category (only applies to td's because Fred case hasn't been added yet)
        if len( height) != 0:


            # make figure before loop
            fig, a0 = plt.subplots(1, 1, figsize=( 8.5, 7) )
            helper_fns.change_font_sizes(small=24, medium=24 )

            # set the figure background (not plot background!) to transparent!!
            fig.patch.set_facecolor('blue')
            fig.patch.set_alpha(0)

            # create a histogram
            plt.sca( a0)

            # sns.displot( x=height, kind='kde')
            sns.set_theme(style="white", palette=None)

            a0.set_xlabel( 'Percentage of CRL Returns at a Given Height')
            a0.set_ylabel( 'Height from Surface (Km)')
            a0.set_ylim( [-.2, 3.75])
            # a0.set_xlim( [0, 200])
            # a0.set_xlim( [-.05, 2])
            a0.grid(True)


            # the line below was probs causing the annoying lack of plot space!!!
            # sns.set_context(rc = {'patch.linewidth': 0.0})
            width = .175
            sns.histplot( y= plot_height, kde=True, binwidth= width, stat="probability", edgecolor='k', linewidth=2, color = 'g')

            a1 = a0.twiny()
            a1.set_xlim( 0, a0.get_xlim()[1] / width)  # similir limits on the y-axis to align the plots
            # a1.xaxis.set_major_formatter( PercentFormatter(1 / width))  # show axis such that 1/binwidth corresponds to 100%
            # a1.set_xlabel( "")

            sns.kdeplot( y= plot_height, bw_adjust= 1, linewidth=2, color = 'b', linestyle='--', ) # bw_adjust=1 -> curve found in histplot

            a0.set_title( "Cloud Height Distribution for " + fig_title_nice[ i], fontsize=16)


            # save the histogram
            os.chdir( "/Users/etmu9498/research/figures/prob-dist-results/intensity-same-plot/")
            if no_eyewalls:
                plt.savefig( fig_title[i] + "-no-eyewalls.png", bbox_inches='tight', dpi=500, transparent=False )
            else:
                plt.savefig( fig_title[i] + ".png", bbox_inches='tight', dpi=500, transparent=False )
            print( fig_title[i] + ' figure created.')

        else:
            print( "The tropical depression case hasn't yet been added!")
            continue
