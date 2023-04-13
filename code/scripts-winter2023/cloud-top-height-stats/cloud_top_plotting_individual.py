# import...
import matplotlib.pyplot as plt
import numpy as np
import os
import xarray as xr
import sys
from scipy.signal import find_peaks
import pandas as pd
# import matplotlib.ticker as ticker # used to manually set colorbar tick number


os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
sys.path.append(  "/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops


# make histogram plots for each individual pass, with the lidar data beneath it!
# savefig = false: save figure to default location and don't enlarge fonts
# savefig = true: save to colloquium slides with high resolution and enlarge fonts!
def plot( tc='all', savefig=False):
    # plot crl data underneath it to highlight how the distributions are created
    metadata = eyewall_metadata.all_metadata()
    lw = 2
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"

    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']

        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( make_plots.load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( make_plots.load_flight_level( crl_data_root + '2022', print_files=False) )

    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2022', print_files=False)]

    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")



    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))

    print( 'data saved to: /Users/etmu9498/research/figures/CRL-all-data-processed/"year"-distributions')

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):

        for filei, fileval in enumerate( filelist[ yeari]):

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

            # make plots for proposed eyewall limits!
            for pairi, pairval in enumerate( eyewall_limits):

                print( "Eyewall limits: " + str( pairval))
                # only do this if there are eyewall limits to look at
                if len( pairval) > 1:
                    lefteyewall = pairval[ 0]
                    righteyewall = pairval[ 1]
                # otherwise, give the eyewalls some arbitrary values and don't trim the eyewalls
                else:
                    lefteyewall = -100
                    righteyewall = -50

                crl_path = crl_data_root + yearval
                os.chdir( crl_path)
                crl_data = xr.open_dataset( fileval)


                ################
                # make the figure!
                fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]}, figsize=(14, 14))
                helper_fns.change_font_sizes( 16, 16)
                time = crl_data.time
                height = crl_data.height
                a0, a1 = ax[0], ax[1]


                ######################
                # plot #1
                # plot power ch1
                # also find and plot cloud top heights!!
                if yearval == '2021':
                    min = -30
                elif yearval == '2022':
                    min = -40

                # preparation
                li = np.argmin( np.abs(crl_data.time.values - lefteyewall))
                ri = np.argmin( np.abs(crl_data.time.values - righteyewall))
                H = crl_data.height
                power = crl_data.P_ch1[ li:ri, :]
                axis = crl_data.time[ li:ri]
                p3_height = crl_data.p3_height[ li:ri]
                cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = min)

                # trim power data for faster plotting / zoomed in on eye
                # find inds 30 mins before li and 30 mins after ri
                shift = .15
                li_new = np.argmin( np.abs(crl_data.time.values - lefteyewall )) # + shift))
                ri_new = np.argmin( np.abs(crl_data.time.values - righteyewall )) # - shift))
                time = time[ li_new:ri_new]
                power = crl_data.P_ch1[ li_new:ri_new]

                print( len( time))
                max = -10

                a1.set_title(str(date[0:2] + '/' + date[2:4] + " Eye Profile"))

                helper_fns.change_font_sizes(16, 16)

                if len( time) > 0:
                    p = a1.pcolormesh( time, height, power.transpose(), vmin = min, vmax = max)
                    a1.plot( cloudtime, cloudheights, c='r', linewidth=lw, label="cloud tops")
                    cbar = plt.colorbar(label="Power Ch. 1 (dBz)", mappable= p, ticks=[-10, -15, -20, -25, -30])
                    # Set the number of ticks
                    # cbar.locator = ticker.MaxNLocator(nbins=5)
                    # cbar.update_ticks()

                a1.legend(loc='upper right')
                a1.set_ylabel("Height (m)", fontsize=16)
                a1.set_facecolor('k')
                a1.set_xlabel( "Time (Hours, UTC)", fontsize=16)

                # plot proposed eyewall limits!
                # a1.axvline( x = lefteyewall, c='w', linewidth=lw)
                # a1.axvline( x = righteyewall, c='w', linewidth=lw)
                a1.set_ylim( [ np.nanmin( height), np.nanmax( height)])
                plt.yticks([1000, 2000, 3000, 4000], fontsize=16)
                plt.xticks(fontsize=16)



                ##################
                # plot #2
                # make a weighted histogram!
                # find the closest left and right indices
                height = cloudheights
                binwidth=25
                smoothwidth=25
                lw=2
                # remove 0 km heights from figure!
                plot_height = height[ np.where( height > 50)[0] ]
                height_bin=np.arange(0, 4510, binwidth)

                # calculate mean height variation within the bin range
                mean_height=[]
                mean_height_count=[]
                mean_height_prob = []

                # do this for every height level determined by the manually inputed bin width
                for newi in height_bin:

                    # find the points that fall within this height bin for this step
                    res=np.where(np.logical_and( plot_height >= newi - binwidth / 2., plot_height <= newi + binwidth / 2. ))

                    # print( plot_height[ res])

                    mean_height.append( np.mean( plot_height[ res]))
                    mean_height_count.append( len( res[0]))
                    # use height, not plot_height, to include heights below 50m in caclulations
                    # this line accounts for clear air fraction when scaling the curves!
                    if len( height) > 0:
                        mean_height_prob.append( len( res[0] ) / len( height))
                    else:
                        # this should only happen if there's no data at this height bin... unlikely,
                        # unless there are no cases at this cat
                        if newi == 0:
                            print( "Divide by 0 attempted :/")
                        mean_height_prob.append( len( res[0] ) / .00001 )

                # smooth data before plotting to eliminate noise
                box_pts = smoothwidth
                box = np.ones(box_pts)/box_pts
                prob_smooth = np.convolve( mean_height_prob, box, mode='same')


                ##########
                ## new code
                # find the color corresponding to the current intensity!
                if date == '0812':
                    col = 'b'
                else:
                    cat = metadata[ yearval]['category'][ date]
                    if cat == 'td':
                        col = 'b'
                    elif cat == 'ts':
                        col = 'k'
                    elif cat == 'wh':
                        col = 'y'
                    elif cat == 'sh':
                        col = 'g'


                # add the probability plot to the total figure! will add nice touches and save later
                a0.plot( prob_smooth * 100, height_bin, color=col, linewidth=lw)

                # add some nice labels, scales, etc
                a0.set_ylabel("Height (m)")
                a0.set_xlabel("Cloud Height Probability (%)")
                a0.set_title( "Cloud Height Probabilities for " + fileval)

                if savefig:
                    # save the total figure!
                    savedir = "/Users/etmu9498/research-private/colloquium"
                    os.chdir( savedir)
                    plt.savefig( fileval[:-3] + '-' + str( pairi) + ".png", dpi=250, bbox_inches='tight')
                else:
                    # save the total figure!
                    savedir = "/Users/etmu9498/research/figures/CRL-all-data-processed/" + yearval + "-distributions"
                    os.chdir( savedir)
                    plt.savefig( fileval[:-3] + '-' + str( pairi) + ".png", dpi=100, bbox_inches='tight')
                print( "Figure " + fileval + ", case " + str( pairi) + " Saved")


# similar to the code above, but don't include the lidar data plots, and plot all
# histograms for a day atop one another
def plot_one_day( tc='all', savefig=False):
    metadata = eyewall_metadata.all_metadata()
    lw = 2
    crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"
    window = 5

    # case 1: do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( make_plots.load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( make_plots.load_flight_level( crl_data_root + '2022', print_files=False) )
    # case 2: do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ make_plots.load_flight_level( crl_data_root + '2022', print_files=False)]
    # case 3: do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")

    # print out the number of files to be saved
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be plotted: " + str( filecount))
    print( 'data saved to: /Users/etmu9498/research/figures/cloud-heights-one-day/"year"')

    old_date = False # date from previous run
    oldyear = False # year from previous run
    # plot colors
    colors = [ 'darkgreen', 'seagreen', 'springgreen', 'mediumturquoise', 'darkcyan']
    # recursively save height data here! for making histograms for a full day
    day_heights = []
    oldname = False

    # make a pandas dataframe to hold number of peaks per distribution per pass info! Will have the following structure:
    # year   |   date   |   pass   |   peaks  |
    # -----------------------------------------
    #  ...   |    ...   |   ....   |    ...   |
    df_peaks = pd.DataFrame( )
    # save values in temporary lists here: will be added to df_peaks once filled
    df_yearlist = []
    df_datelist = []
    df_passlist = []
    df_peaklist = []
    df_peakcount_strong = []
    df_peakcount_weak = []
    df_intensity = []
    df_category = []


    # do this for all selected datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, fileval in enumerate( filelist[ yeari]):
            ######
            ## new code: grab the limits for this case
            ######
            date = fileval[7:11]

            # check if this date exists... if not, give it some empty eyewall limits!
            # also check for fred am and pm cases!!
            if date == '0812':
                if fileval[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                elif fileval[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
            else:
                eyewall_limits = [ ()]

            # make plots for proposed eyewall limits!
            for pairi, pairval in enumerate( eyewall_limits):
                print( "Eyewall limits: " + str( pairval))
                # only do this if there are eyewall limits to look at
                if len( pairval) > 1:
                    lefteyewall = pairval[ 0]
                    righteyewall = pairval[ 1]
                # otherwise, give the eyewalls some arbitrary values and don't trim the eyewalls
                else:
                    lefteyewall = -100
                    righteyewall = -50
                crl_path = crl_data_root + yearval
                os.chdir( crl_path)
                crl_data = xr.open_dataset( fileval)

                time = crl_data.time
                height = crl_data.height
                if yearval == '2021':
                    min = -30
                elif yearval == '2022':
                    min = -40
                # preparation: find cloud heights for this eye pass
                li = np.argmin( np.abs(crl_data.time.values - lefteyewall))
                ri = np.argmin( np.abs(crl_data.time.values - righteyewall))
                H = crl_data.height
                power = crl_data.P_ch1[ li:ri, :]
                axis = crl_data.time[ li:ri]
                p3_height = crl_data.p3_height[ li:ri]
                cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = min)

                # first new date case: no need to save the old figure, just make a new one
                if not old_date:
                    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
                    if savefig:
                        helper_fns.change_font_sizes( 18, 18)
                    else:
                        helper_fns.change_font_sizes( 12, 12)

                    # add some nice labels, scales, etc
                    ax.set_ylabel("Height (m)")
                    ax.set_xlabel("Cloud Height Probability (%)")
                    ax.set_title( "Cloud Height Probabilities for " + fileval)

                # if it's a new date (but not the first), save the old figure and make a new one!
                elif old_date != date:
                    count = 0 # reset the color counter
                    # saving the multiple distribution figure for the previous date

                    # add scatter points again just to add labels! only once
                    # add scatter points symbolyzing tall / strong peak locations!
                    ax.scatter( prob_smooth[ peakind], height_bin[ peakind], c='k', marker='*', s=50, zorder=5, label='strong')
                    # add scatter points symbolyzing weaker / smaller peak locations!
                    ax.scatter( prob_smooth[ peakindweak], height_bin[ peakindweak], c='r', marker='*', s=50, zorder=6, label='weak')

                    ax.legend( loc='upper right')
                    savedir = "/Users/etmu9498/research/figures/cloud-heights-one-day/" + oldyear
                    os.chdir( savedir)
                    fig.savefig( old_date + ".png", dpi=100, bbox_inches='tight')
                    print( "Figure " + old_date + " Saved")

                    # add a new plot with histogram data averaged from all passes for one date!
                    # get the smoothed probability function
                    prob_smooth, height_bin = weighting_fn( np.array( day_heights))



                    # *** # 
                    # new code: smooth again to get rid of big height bumps?! maybe don't plot like this, but 
                    # useful for finding bimodal, etc peaks
                    prob_smooth = np.array( pd.Series( prob_smooth).rolling(window=window, min_periods=1, center=True).mean()) 

                    # create a new many subplot figure and add histogram data
                    fig2, ax2 = plt.subplots(1, 1, figsize=(9, 6))
                    helper_fns.change_font_sizes( 12, 12)

                    # add some nice labels, scales, etc
                    ax2.set_ylabel("Height (m)")
                    ax2.set_xlabel("Cloud Height Probability (%)")
                    ax2.set_title( "Cloud Height Probabilities for " + oldname)
                    # find the color corresponding to the current intensity for plotting!
                    if old_date == '0812':
                        col = 'b'
                    else:
                        cat = metadata[ oldyear]['category'][ old_date]
                        if cat == 'td':
                            col = 'b'
                        elif cat == 'ts':
                            col = 'k'
                        elif cat == 'wh':
                            col = 'y'
                        elif cat == 'sh':
                            col = 'g'
                    # add the total probability plot to the figure!
                    ax2.plot( prob_smooth, height_bin, color=col, linewidth=lw)

                    # find local peaks for this probability
                    peakind, peakindweak = find_local_peak( prob_smooth)

                    # add scatter points symbolyzing tall / strong peak locations!
                    ax2.scatter( prob_smooth[ peakind], height_bin[ peakind], c='k', marker='*', s=50, zorder=5, label='strong')
                    # add scatter points symbolyzing weaker / smaller peak locations!
                    ax2.scatter( prob_smooth[ peakindweak], height_bin[ peakindweak], c='r', marker='*', s=50, zorder=6, label='weak')

                    savedir = "/Users/etmu9498/research/figures/cloud-heights-one-day/" + oldyear
                    os.chdir( savedir)
                    fig2.savefig( old_date + "-total.png", dpi=100, bbox_inches='tight')
                    print( "Figure " + old_date + "-total Saved")
                    # reset the day_heights after creating this day's histogram!
                    day_heights = []

                    # create a new many subplot figure and add histogram data
                    fig, ax = plt.subplots(1, 1, figsize=(9, 6))
                    helper_fns.change_font_sizes( 12, 12)
                    # add some nice labels, scales, etc
                    ax.set_ylabel("Height (m)")
                    ax.set_xlabel("Cloud Height Probability (%)")
                    ax.set_title( "Cloud Height Probabilities for " + fileval)

                # use the helper function to weight and smooth the crl cloud heights!
                prob_smooth, height_bin = weighting_fn( cloudheights)

                # print( prob_smooth)
                # print( type( prob_smooth))

                # new code from above
                prob_smooth = np.array( pd.Series( prob_smooth).rolling(window=window, min_periods=1, center=True).mean())

                # add the probability plot to the figure!
                ax.plot( prob_smooth, height_bin, color=colors[ pairi], linewidth=lw, label='case ' + str( pairi))

                # find local peaks for this probability
                peakind, peakindweak = find_local_peak( prob_smooth)

                # *** # 
                # add the metadata and peaks found to the dataframe!
                df_yearlist.append( yearval)
                df_datelist.append( date)
                df_passlist.append( pairi)
                df_peaklist.append( height_bin[peakind].tolist())
                df_peakcount_strong.append( len( height_bin[peakind]))
                df_peakcount_weak.append( len( height_bin[peakindweak]))
                
                if date == '0812':
                    if fileval[11:13] == "H1":
                        df_intensity.append( metadata[yearval]['intensity']['0812am'])
                        df_category.append( metadata[yearval]['category']['0812am'])
                    elif fileval[11:13] == "H2":
                        df_intensity.append( metadata[yearval]['intensity']['0812pm'])
                        df_category.append( metadata[yearval]['category']['0812pm'])
                else:
                    df_intensity.append( metadata[yearval]['intensity'][date])
                    df_category.append( metadata[yearval]['category'][date])

                # add scatter points symbolyzing tall / strong peak locations!
                ax.scatter( prob_smooth[ peakind], height_bin[ peakind], c='k', marker='*', s=50, zorder=5)
                # add scatter points symbolyzing weaker / smaller peak locations!
                ax.scatter( prob_smooth[ peakindweak], height_bin[ peakindweak], c='r', marker='*', s=50, zorder=6)

                # update the old date / year counter, add one to count to iterate colors
                old_date = date
                oldyear = yearval
                oldname = fileval

                # add the heights to the list for this date!
                day_heights += cloudheights.tolist()

    # made it outside of year / day loops: save the last remaining averaged and many height plots now
    # saving the multiple distribution figure for the previous date
    ax.legend( loc='upper right')
    savedir = "/Users/etmu9498/research/figures/cloud-heights-one-day/" + oldyear
    os.chdir( savedir)
    fig.savefig( old_date + ".png", dpi=100, bbox_inches='tight')
    print( "Figure " + old_date + " Saved")

    # add a new plot with histogram data averaged from all passes for one date!
    prob_smooth, height_bin = weighting_fn( np.array( day_heights))
    # new code from above
    prob_smooth = np.array( pd.Series( prob_smooth).rolling(window=window, min_periods=1, center=True).mean())

    # create a new many subplot figure and add histogram data
    fig2, ax2 = plt.subplots(1, 1, figsize=(9, 6))
    helper_fns.change_font_sizes( 12, 12)
    # add some nice labels, scales, etc
    ax2.set_ylabel("Height (m)")
    ax2.set_xlabel("Cloud Height Probability (%)")
    ax2.set_title( "Cloud Height Probabilities for " + oldname)
    # find the color corresponding to the current intensity for plotting!
    if old_date == '0812':
        col = 'b'
    else:
        cat = metadata[ oldyear]['category'][ old_date]
        if cat == 'td':
            col = 'b'
        elif cat == 'ts':
            col = 'k'
        elif cat == 'wh':
            col = 'y'
        elif cat == 'sh':
            col = 'g'
    # add the total probability plot to the figure!
    ax2.plot( prob_smooth, height_bin, color=col, linewidth=lw)

    # find local peaks for this probability
    peakind, peakindweak = find_local_peak( prob_smooth)

    # add scatter points symbolyzing tall / strong peak locations!
    ax.scatter( prob_smooth[ peakind], height_bin[ peakind], c='k', marker='*', s=50, zorder=5, label='strong')
    # add scatter points symbolyzing weaker / smaller peak locations!
    ax.scatter( prob_smooth[ peakindweak], height_bin[ peakindweak], c='r', marker='*', s=50, zorder=6, label='weak')

    savedir = "/Users/etmu9498/research/figures/cloud-heights-one-day/" + oldyear
    os.chdir( savedir)
    fig2.savefig( old_date + "-total.png", dpi=100, bbox_inches='tight')
    print( "Figure " + old_date + "-total Saved")


    # finally, assemble and return the distribution peak dataframe!


    df_peaks['year'] = df_yearlist
    df_peaks['date'] = df_datelist
    df_peaks['pass'] = df_passlist
    df_peaks['peak heights (m)'] = df_peaklist
    df_peaks['peak counts strong'] = df_peakcount_strong
    df_peaks['peak counts weak'] = df_peakcount_weak
    df_peaks['intensity (kt)'] = df_intensity
    df_peaks['tc category'] = df_category

    return df_peaks


# helper function for finding distribution peaks! 
# change this function to update peaks for the 4 or 5 times this is called above
def find_local_peak( prob_smooth):
    # max height of the distribution
    maxh = np.nanmax( prob_smooth)
    # no height requirements, but the peaks must be prominent enough!
    peakind = find_peaks( x= prob_smooth, prominence= maxh / 6, distance=20)[0].tolist()

    print("tall: " + str( peakind))

    # remove tall peaks from this list!
    peakindweak = find_peaks( x= prob_smooth, prominence= maxh / 15, distance=10)[0].tolist()
    peakindweak = list( set(peakind).symmetric_difference( set(peakindweak)))

    # get rid of the strong peaks to just look at the weak peaks!
    print("Smaller: " + str( peakindweak))



    return peakind, peakindweak



# helper function for making a weighted histogram!
def weighting_fn( height):
    binwidth=25
    smoothwidth=25
    lw=1.5
    # remove 0 km heights from figure!
    plot_height = height[ np.where( height > 50)[0] ]
    height_bin=np.arange(0, 4510, binwidth)

    # calculate mean height variation within the bin range
    mean_height=[]
    mean_height_count=[]
    mean_height_prob = []
    # do this for every height level determined by the manually inputed bin width
    for newi in height_bin:
        # find the points that fall within this height bin for this step
        res=np.where(np.logical_and( plot_height >= newi - binwidth / 2., plot_height <= newi + binwidth / 2. ))
        mean_height.append( np.mean( plot_height[ res]))
        mean_height_count.append( len( res[0]))
        # use height, not plot_height, to include heights below 50m in caclulations
        # this line accounts for clear air fraction when scaling the curves!
        if len( height) > 0:
            mean_height_prob.append( len( res[0] ) / len( height))
        else:
            # this should only happen if there's no data at this height bin... unlikely,
            # unless there are no cases at this cat
            if newi == 0:
                print( "Divide by 0 attempted :/")
            mean_height_prob.append( len( res[0] ) / .00001 )

    # smooth data before plotting to eliminate noise
    box_pts = smoothwidth
    box = np.ones(box_pts)/box_pts
    prob_smooth = np.convolve( mean_height_prob, box, mode='same')

    return prob_smooth * 100, height_bin
