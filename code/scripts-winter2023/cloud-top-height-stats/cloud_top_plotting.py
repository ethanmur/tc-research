# import...
import matplotlib.pyplot as plt
import numpy as np
import os

def intensity_bins( year, td_heights, ts_heights, wh_heights, sh_heights, binwidth=25, smoothwidth=25, lw=2):
    # create histograms for every TC intensity
    # put things in lists for easier looping
    fig_title = [ 'tropical-depressions', 'tropical-storms', 'weak-hurricanes', 'strong-hurricanes']
    # nicely formatted titles for plots, as opposed to saving
    fig_title_nice = [ 'Tropical Depressions', 'Tropical Storms', 'Weak Hurricanes', 'Strong Hurricanes']
    colors = [ 'b', 'k', 'y', 'g']
    heights = [ td_heights, ts_heights, wh_heights, sh_heights]
    savedir = "/Users/etmu9498/research/figures/new-intensity-dists/"

    # make figure before loop
    fig, a0 = plt.subplots(1, 1, figsize=( 7, 7) )
    smallfont=16
    mediumfont=16
    largefont=20
    # set the figure background (not plot background!) to transparent!!
    # fig.patch.set_facecolor('blue')
    # fig.patch.set_alpha(0)

    # create a histogram
    plt.sca( a0)
    # loop through each case
    for i in range( 4):
        height = np.array( heights[ i])
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

        # remove 0% probability lines (mostly at high alts) from histogram plots!
        prob_smooth_trim = prob_smooth[ np.where( prob_smooth > 0.00001)[0]]
        height_trim = height_bin[ np.where( prob_smooth > 0.00001)[0]]

        # add the probability plot to the total figure! will add nice touches and save later
        a0.plot( prob_smooth_trim * 100, height_trim, color=colors[i], label=fig_title_nice[i], linewidth=lw)


        # use the highest smoothed value to effectively scale the x axis!!
        max_prob = np.nanmax( prob_smooth)
        # print( 'max prob: ' + str( max_prob))


        # make and save separate plots for each intensity! in the same output folder
        fig1, a1 = plt.subplots(1, 1, figsize=( 7, 7) )
        # set the figure background (not plot background!) to transparent!!
        fig1.patch.set_facecolor('blue')
        fig1.patch.set_alpha(0)
        plt.sca( a1)
        a1.plot( np.array( mean_height_prob) * 100, height_bin, color=colors[i], label='Raw Data: Bins = ' + str( binwidth * 1000) + " m", linewidth=lw/2, alpha=.7)
        a1.plot( prob_smooth * 100, height_bin, color=colors[i], label='Smoothed Data', linewidth=lw)

        # a1.set_xlabel('Cloud Height Probability (Sums to 1)', fontsize=mediumfont)
        a1.set_xlabel('Cloud Height Probability (%)', fontsize=mediumfont)

        # sns.set_theme(style="white", palette=None)
        a1.set_ylabel( 'Height from Surface (Km)', fontsize=mediumfont)
        # a1.set_ylim( [-.2, 4.5])
        a1.set_ylim( [-200, 4500])

        # a0.set_xlim( [ 0, 1.05])
        # a1.set_xlim( [ 0, .02])
        a1.set_xlim( [ 0, max_prob * 2 * 100])


        a1.grid(True)
        a1.set_title( "Cloud Height Distributions for " + fig_title_nice[i], fontsize=largefont)
        a1.tick_params(axis='both', which='major', labelsize=smallfont)
        a1.legend(fontsize=smallfont, loc="upper right")

        os.chdir( savedir)
        plt.savefig( year + "-" + fig_title[i] + '-' + str( int( binwidth*1000)) + "m.png", bbox_inches='tight', dpi=500, transparent=False )
        plt.sca( a0)

    # save the total figure!
    a0.set_ylabel( 'Height from Surface (Km)', fontsize=mediumfont)
    a0.set_xlabel('Cloud Height Probability (%)', fontsize=mediumfont)
    a0.set_ylim( [-200, 4600])
    a0.grid(True)
    a0.set_title( "Cloud Height Distributions by TC Intensity", fontsize=largefont)
    a0.tick_params(axis='both', which='major', labelsize=smallfont)
    plt.legend(fontsize=smallfont, loc="upper right")
    plt.savefig( year + "-all-dists-" + str( int( binwidth*1000)) + "m.png", bbox_inches='tight', dpi=500, transparent=False )


# make a simple plot of eye cloud top heights!
def eye_cloud_tops( heights, time, axis, H, power, title):

    plt.figure( figsize=(10, 5))
    plt.title( title)
    min = -30
    max = -10
    plt.pcolormesh( axis, H, power.transpose(), vmin = min, vmax = max)
    plt.plot( time, heights, c='r')
    plt.xlabel( 'Time (UTC, Hours)')
    plt.ylabel( "Height (km)")
    plt.ylim( [0, np.nanmax( H)])
    ax = plt.gca()
    ax.set_facecolor('k')
