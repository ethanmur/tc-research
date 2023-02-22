import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import xarray as xr
import math
import seaborn as sns
import scipy

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import make_plots_new_heights
import helper_fns
import cloud_height

os.chdir(  "/Users/etmu9498/research/code/scripts/plotting")
import simple_flight_level_plot


# plot and save a panel of height distributions vs rmw, crl, and in situ data for each
# TC eye pass for the requested tc
# valid test types:
# == 'tc_plots': calc stats and plot crl eye passes + stats for each case. Takes a while
# == 'intensity': calc stats for each intensity category and make height plots for
#    different intensities. Pretty fast.
# == 'stats': only calc stats, don't make plots. Really fast
# == 'return': return arrays holding height values for whole, inner, and outer cases
#    for analysis outside of this script!
def plot_all( tc='all', print_stats=True, test_type='tc_plots', binwidth=.025, smoothwidth=20, plot_raw=False):

    # code taken from cloud_height_pdfs_in_situ.pdf_all_eyes()
    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['fred', 'grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    # initializing lists to hold TC cloud height vs intensity data for test_type='intensity'!
    td_heights, td_inner, td_outer, td_cases = [], [], [], 0
    ts_heights, ts_inner, ts_outer, ts_cases = [], [], [], 0
    wh_heights, wh_inner, wh_outer, wh_cases = [], [], [], 0
    sh_heights, sh_inner, sh_outer, sh_cases = [], [], [], 0

    # lists for test_type='tc_plots'
    # an empty list that will hold all the height datasets from within the eye
    # this is done in the same way as the wrapper method found in cloud_height_pdfs_in_situ.py
    all_heights = []

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

        # setup
        tcname = metadata['tc_name']
        # do this for every eye pass for this tc
        for dataset in range( len( metadata[ 'dates'] )):

            # default, old test: calc stats and make figures for individual TCs
            if test_type == 'tc_plots':
                # call the helper function to find a new rmw axis and create a panelled plot
                newH = plot_one( tcname, dataset, print_stats)
                savedir = "/Users/etmu9498/research/figures/rmw-inner-vs-outer"
                # add new calculated heights to previously found heights for this case
                all_heights = all_heights + np.ndarray.tolist( newH)

                # save the newly created figure
                os.chdir( savedir)
                plt.savefig( metadata['tc_name'].casefold() + "-" + str( dataset+1) + ".png", bbox_inches='tight', dpi=200 )
                print( "Plot " + str( dataset + 1) + " saved\n" )

            # find stats based on tc intensity
            elif test_type == 'intensity' or test_type == 'return':
                # call a shorter helper function to just find the cloud heights for this case: no plotting yet
                H, lowrmw, highrmw = findH( tcname, dataset)

                # save height values from this run
                cat = metadata['intensity_cat'][dataset] # tc intensity category
                # figure out where to put the height data depending on the tc intensity!
                if cat == 'td':
                    td_heights += np.ndarray.tolist( H)
                    td_inner += lowrmw
                    td_outer += highrmw
                    # print( 'td')
                    td_cases += 1
                elif cat == 'ts':
                    # print( 'ts')
                    ts_heights += np.ndarray.tolist( H)
                    ts_inner += lowrmw
                    ts_outer += highrmw
                    ts_cases += 1
                elif cat == 'wh':
                    # print( 'wh')
                    wh_heights += np.ndarray.tolist( H)
                    wh_inner += lowrmw
                    wh_outer += highrmw
                    wh_cases += 1
                elif cat == 'sh':
                    # print( 'sh')
                    sh_heights += np.ndarray.tolist( H)
                    sh_inner += lowrmw
                    sh_outer += highrmw
                    sh_cases += 1


            elif test_type == 'stats':
                # find heights for this run and calc cloud top height stats
                H, lowrmw, highrmw = findH( tcname, dataset)
                H = np.array( H)
                calc_stats( H)

                # save current eye pass stats in list to eventually look at all stats
                all_heights = all_heights + np.ndarray.tolist( H)
                print( "Stats for eye pass " + str( dataset + 1) + " calculated\n" )

            else:
                print("Please choose a valid entry for test_type: 'tc_plots', 'intensity', 'return', or 'stats'.")
                return


        # this snippet is outside the previous inner loop. should be completed after all the other steps!
        # This should be run when we get to the end of one tc; should be done for henri, grace, etc in order
        if test_type == 'tc_plots' or test_type == 'stats':
            # once done cycling through every case for one tc, print out total stats from that case
            print( "Overall Statistics for All Eye Passes, TC " + metadata['tc_name'] + ":")
            all_heights = np.array( all_heights)
            calc_stats( all_heights)

            # make all heights an empty python list again! to prevent duplicatations
            all_heights = []


    # this check is even further outside the previous statement! Create intensity plots
    # once height finder scripts have been run for all tc eye passes
    if test_type == 'intensity':
        # turn assorted variables into lists for easier data transfer
        cases = [ td_cases, ts_cases, wh_cases, sh_cases]
        heights = [ td_heights, ts_heights, wh_heights, sh_heights]
        inner = [ td_inner, ts_inner, wh_inner, sh_inner]
        outer = [ td_outer, ts_outer, wh_outer, sh_outer]

        # calc stats for each intensity category!
        category = ["TDs:", "TSs:", "WHs:", "SHs:"]
        for i in range( 4):
            print( "\nAll height statistics for " + category[ i])
            calc_stats( np.array( heights[ i] ), c=.75)
            print('\n')
            calc_stats( np.array( heights[ i] ), c=.95)
            print( '\n')
            calc_stats( np.array( heights[ i] ), c=.99)

        # make intensity figures with this step
        plot_intensity( cases, heights, inner, outer, binwidth, smoothwidth, plot_raw)
        # plot_all_inner_outer( cases, heights, inner, outer)

    elif test_type == 'return':
        heights = [ td_heights, ts_heights, wh_heights, sh_heights]
        inner = [ td_inner, ts_inner, wh_inner, sh_inner]
        outer = [ td_outer, ts_outer, wh_outer, sh_outer]

        return heights, inner, outer


# helper function used to create plots of dists of cloud heights by distance from TC center: total, inner, and outer.
# the total plot should be the same as the one created for the esss poster! It's mostly
# there to check work.
def plot_all_inner_outer( cases, heights, inner, outer):

    int_list = ["Tropical Depressions", "Tropical Storms", "Weak Hurricanes", "Strong Hurricanes"]
    color_list = [ 'b', 'k', 'y', 'g']
    padding = 100
    plt.figure(figsize=(24, 7))
    helper_fns.change_font_sizes(small=20, medium=20 )
    width = .175 # .3  1.0
    lw = 2.0

    # add each intensity distribution to the correct subplot
    for i in range( 4):

        # choose the right xaxis groups for this intensity case, much like before.
        # sorry for the confusing syntax :(
        height = np.array( heights[ i])
        inneri = np.array( inner[ i])
        outeri = np.array( outer[ i])

        # remove 0 km heights from figure! Precisely, they need to be above 50 m
        height = height[ np.where( height > .05)[0] ]
        inneri = inneri[ np.where( inneri > .05)[0] ]
        outeri = outeri[ np.where( outeri > .05)[0] ]

        # normal early case:
        if i < 3:
            # add histograms!
            plt.subplot(131)
            sns.kdeplot( y= height, bw_adjust= width, linewidth=lw, color = color_list[ i], cut=0)
            plt.subplot(132)
            sns.kdeplot( y= inneri, bw_adjust= width, linewidth=lw, color = color_list[ i], cut=0)
            plt.subplot(133)
            sns.kdeplot( y= outeri, bw_adjust= width, linewidth=lw, color = color_list[ i], label=int_list[ i], cut=0)

        # last case: add labels, legend titles, etc
        else:
            # add histograms!
            plt.subplot(131)
            sns.kdeplot( y= height, bw_adjust= width, linewidth=lw, color = color_list[ i], cut=0)
            plt.xlabel( 'Cloud Height Density')
            plt.ylabel( 'Height from Surface (Km)')
            plt.title( "Height Distribution: Whole Eye")
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, 1.5])
            plt.grid(True)

            plt.subplot(132)
            sns.kdeplot( y= inneri, bw_adjust= width, linewidth=lw, color = color_list[ i], cut=0)
            plt.xlabel( 'Cloud Height Density')
            plt.ylabel( 'Height from Surface (Km)')
            plt.title( "Height Distribution: Inner Eye")
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, 1.5])
            plt.grid(True)

            plt.subplot(133)
            sns.kdeplot( y= outeri, bw_adjust= width, linewidth=lw, color = color_list[ i], label=int_list[ i], cut=0)
            plt.xlabel( 'Cloud Height Density')
            plt.ylabel( 'Height from Surface (Km)')
            plt.title( "Height Distribution: Outer Eye")
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, 1.5])
            plt.grid(True)

            plt.legend(loc='lower right', fontsize = 15)

    # save figure
    savedir = "/Users/etmu9498/research/figures/rmw-inner-vs-outer/intensity-new"
    os.chdir( savedir)
    save_list = ["tropical-depressions", "tropical-storms", "weak-hurricanes", "strong-hurricanes"]
    # plt.savefig( save_list[i] + ".png", bbox_inches='tight', dpi=200 )





# helper function to create updated intensity plots! There's a lot more detail in these plots
def plot_intensity( cases, heights, inner, outer, binwidth=.1, smoothwidth=20, plot_raw=False):
    # loop through each intensity case

    for i in range( 4):
        height = np.array( heights[ i])
        inneri = np.array( inner[ i])
        outeri = np.array( outer[ i])

        int_list = ["Tropical Depressions", "Tropical Storms", "Weak Hurricanes", "Strong Hurricanes"]


        # remove 0 km heights from figure! Precisely, they need to be above 50 m
        height = height[ np.where( height > .05)[0] ]
        inneri = inneri[ np.where( inneri > .05)[0] ]
        outeri = outeri[ np.where( outeri > .05)[0] ]


        # check if there are no cases for this category (not an issue after adding fred td passes)
        if len( height) != 0:

            padding = 100
            plt.figure(figsize=(8, 8))
            helper_fns.change_font_sizes(small=14, medium=14 )

            width = binwidth
            lw = 2.0

            # add histograms!
            # smoothing code taken from cloud_height_pdfs_all_one_figure.py
            height_bin=np.arange(0, 4.51, width)
            inner_height_prob = []
            outer_height_prob = []


            for newi in height_bin:
                inner_res=np.where(np.logical_and( inneri >= newi - width / 2., inneri <= newi + width / 2. ))
                outer_res=np.where(np.logical_and( outeri >= newi - width / 2., outeri <= newi + width / 2. ))

                inner_height_prob.append( len( inner_res[0] ) / len( height))
                outer_height_prob.append( len( outer_res[0] ) / len( height))


            # smooth data before plotting to eliminate noise
            box_pts = smoothwidth
            box = np.ones(box_pts)/box_pts
            inner_smooth = np.convolve( inner_height_prob, box, mode='same')
            outer_smooth = np.convolve( outer_height_prob, box, mode='same')



            plt.xlabel('Cloud Height Probability (%)')


            # sns.kdeplot( y= height, bw_adjust= width, linewidth=lw, color = 'k', label='all clouds', cut=0) #  clip=[0, 3.6]) # bw_adjust=1 -> curve found in histplot
            # plt.xlabel( 'Cloud Height Density (Area Under Curve = 1)')
            plt.ylabel( 'Height from Surface (Km), Flat Lines = Mean')

            plt.title( "Cloud Height Distribution: " + int_list[i])
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, 1.5])

            plt.grid(True)

            plt.plot( outer_smooth * 100, height_bin, color='g', label='outer clouds', linewidth=lw)
            plt.plot( inner_smooth * 100, height_bin, color='b', label='inner clouds', linewidth=lw)

            prob_list = [ np.nanmax( outer_smooth), np.nanmax(inner_smooth) ]
            max_prob = np.nanmax( np.array( prob_list))
            plt.xlim( [ 0, max_prob * 2 * 100])


            # sns.kdeplot( y= inneri, bw_adjust= width, linewidth=lw, color = 'b', label='inner clouds', cut=0) #  clip=[0, 3.6]) # bw_adjust=1 -> curve found in histplot
            # sns.kdeplot( y= outeri, bw_adjust= width, linewidth=lw, color = 'g', label='outer clouds', cut=0) #  clip=[0, 3.6]) # bw_adjust=1 -> curve found in histplot
            plt.legend(loc='lower right')

            # plot mean heights for inner and outer clouds as a horizontal line
            # plt.axhline( y=np.mean( inneri), c='b', linewidth=lw)
            # plt.axhline( y=np.mean( outeri), c='g', linewidth=lw)
            # plt.legend(loc='lower right')

            # option to plot unsmoothed data behind original curves
            save_list = ["tropical-depressions", "tropical-storms", "weak-hurricanes", "strong-hurricanes"]
            if plot_raw:
                plt.plot( np.array( outer_height_prob) * 100, height_bin, color='g', label='outer clouds', alpha = .6, linewidth=lw/2)
                plt.plot( np.array(inner_height_prob) * 100, height_bin, color='b', label='inner clouds', alpha = .6, linewidth=lw/2)
                savename = save_list[i] + "-" + str( int( binwidth * 1000)) + "m-raw.png"
            else:
                savename = save_list[i] + "-" + str( int( binwidth * 1000)) + "m.png"

            # save figure
            savedir = "/Users/etmu9498/research/figures/rmw-inner-vs-outer/intensity-new"
            os.chdir( savedir)
            plt.savefig( savename, bbox_inches='tight', dpi=200 )






# helper function used to create total, inner, and outer cloud height plots for each
# tc intensity!
def plot_intensity_old( cases, heights, inner, outer):
    # loop through each intensity case

    for i in range( 4):
        height = np.array( heights[ i])
        inneri = np.array( inner[ i])
        outeri = np.array( outer[ i])

        int_list = ["Tropical Depressions", "Tropical Storms", "Weak Hurricanes", "Strong Hurricanes"]

        #######################
        # maybe uncomment this section!
        #######################

        # remove 0 km heights from figure!
        height = height[ np.where( height > .05)[0] ]
        inneri = inneri[ np.where( inneri > .05)[0] ]
        outeri = outeri[ np.where( outeri > .05)[0] ]


        # check if there are no cases for this category (not an issue after adding fred td passes)
        if len( height) != 0:
            padding = 100
            plt.figure(figsize=(20, 8))
            helper_fns.change_font_sizes(small=20, medium=20 )

            # add histograms!
            plt.subplot(131)
            sns.histplot( y= height, kde=True, binwidth=.175, stat='probability', edgecolor='k', linewidth=2, color = 'g')
            plt.xlabel( 'Cloud Height Probability (Sums to 1)')
            plt.ylabel( 'Height from Surface (Km)')
            plt.title( "Height Distribution: Whole Eye")
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, .5])
            plt.grid(True)

            title = ( "Cloud Top Height Distributions for " + int_list[ i] )
            axis0 = plt.gca()
            axis0.text( .4, 4.5, title, fontsize=20)

            plt.subplot(132)
            sns.histplot( y= inneri, kde=True, binwidth=.175, stat='probability', edgecolor='k', linewidth=2, color = 'g')
            plt.xlabel( 'Cloud Height Probability (Sums to 1)')
            plt.ylabel( 'Height from Surface (Km)')
            plt.title( "Height Distribution: Inner Eye")
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, .5])
            plt.grid(True)

            plt.subplot(133)
            sns.histplot( y= outeri, kde=True, binwidth=.175, stat='probability', edgecolor='k', linewidth=2, color = 'g')
            plt.xlabel( 'Cloud Height Probability')
            plt.ylabel( 'Height from Surface (Km)')
            plt.title( "Height Distribution: Outer Eye")
            plt.ylim( [-.2, 3.75])
            plt.xlim( [0, .5])
            plt.grid(True)

            # save figure
            savedir = "/Users/etmu9498/research/figures/rmw-inner-vs-outer/intensity"
            os.chdir( savedir)
            save_list = ["tropical-depressions", "tropical-storms", "weak-hurricanes", "strong-hurricanes"]
            plt.savefig( save_list[i] + ".png", bbox_inches='tight', dpi=200 )


            # print out statistics for this case
            # the helper function wasn't used bc I can't pass it an rmw axis: the height data is locally saved :(
            print("Statistics for " + int_list[ i] + ":")
            print( "Number of data points:  " + str( len( height)))
            print( "Height value range:     " + str( height.min()) + " km to " + str( height.max()) + " km")
            print( "Height value median:    " + str( np.median( height)) + " km")
            print( "Standard deviation:     " + str( np.std( height)) )
            print( "Height value mean:      " + str( np.mean( height)) + " km")
            print( "Inner RMW mean H:       " + str( np.mean( inneri)))
            print( "Outer RMW mean H:       " + str( np.mean( outeri)))




# similar to the function below, but it only finds the current height distribution;
# it doesn't calc any stats yet!
# this function also finds and returns inner and outer rmw heights!
def findH( tcname, dataset):

    metadata = tc_metadata.all_data( tcname)
    # load crl data
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    crl_path =  "/Users/etmu9498/research/data/crl-new-matrices"
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # find cutoffs for cloud height calculation
    lim0 = metadata['eyewall_dists_no_eyewalls'][dataset][0]
    lim1 = metadata['eyewall_dists_no_eyewalls'][dataset][1]
    ind0 = np.argmin( (np.abs(crl_data.in_situ_distance - lim0 )).values, axis=0)
    ind1 = np.argmin( (np.abs(crl_data.in_situ_distance - lim1 )).values, axis=0)

    # old code: works the same, but just deprecated and printing out an annoying warning
    # ind0 = (np.abs(crl_data.in_situ_distance - lim0 )).argmin().values
    # ind1 = (np.abs(crl_data.in_situ_distance - lim1 )).argmin().values

    # find cloud heights for current eye
    H, xaxis_value = cloud_height.find_cloud_heights( crl_name, -30, ind0, ind1, xaxis ='in-situ-dist', crl_path=crl_path, new_heights=True)

    ##########################
    ## simpler and probably better way of doing things!! just initialize both sides of the rmw at once
    ##########################
    dist = crl_data.in_situ_distance
    rmwaxis = np.linspace( -1, 1, num=len( dist[ ind0:ind1]))

    # break the distribution down by rmw distance!
    lowrmw = []
    highrmw = []
    # cycle through every case
    for i in range( len( H)):
        # high rmw case: add cloud heights
        if rmwaxis[ i] > .5 or rmwaxis[ i] < -.5:
            highrmw.append( H[ i])
        # low rmw case
        else:
            lowrmw.append( H[ i])
    return H, lowrmw, highrmw



# this helper function actually creates the panelled plot mentioned above. It finds
# a new rmw axis by defining the tc center as the point where there are equal amounts of
# data on either side of the eye, and calls that 0 km. It then finds statistics within
# and outside of rmw = .5 (inner vs outer eye cloud heights)
def plot_one( tcname, dataset, print_stats):

    metadata = tc_metadata.all_data( tcname)
    # load crl data
    tdr_name, crl_name = tc_metadata.choose_new_data( tcname, dataset)
    crl_path =  "/Users/etmu9498/research/data/crl-new-matrices"
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # load flight level data
    flight_name = tc_metadata.choose_new_in_situ_name( tcname, dataset)
    flight_path = "/Users/etmu9498/research/data/in-situ-new"
    os.chdir( flight_path)

    # find the center point of the eye: equal number of data points to the left as to the right
    lim0 = metadata['eyewall_dists_no_eyewalls'][dataset][0]
    lim1 = metadata['eyewall_dists_no_eyewalls'][dataset][1]
    dist = crl_data.in_situ_distance
    # the indices closest to lim0 and lim1
    ind0 = (np.abs(crl_data.in_situ_distance - lim0 )).argmin().values
    ind1 = (np.abs(crl_data.in_situ_distance - lim1 )).argmin().values


    # trim down distance axis to look at the eye
    trimdist = dist[ ind0: ind1]

    # find the eye center
    # even case
    if len( trimdist.values) % 2 == 0:
        small_ind = len( trimdist.values) / 2
    # odd length of xaxis case
    else:
        # math.ceil rounds the division up to get the center int!
        small_ind = math.ceil( len( trimdist.values) / 2)

    # add all the indices before the first eyewall to get a properly spaced tc center
    # the values below are used to plot the tc center below
    startdist = dist[ np.where( dist.values <= lim0)]
    center_ind = int( len( startdist.values) + small_ind )
    center_dist = dist[ center_ind]

    ##########################
    ## simpler and probably better way of doing things!! just initialize both sides of the rmw at once
    ##########################
    rmwaxis = np.linspace( -1, 1, num=len( dist[ ind0:ind1]))

    # find inds closest to rmw = +-.5
    # add distance to ind = 0 for eyewall like above!
    ind_minus5 = (np.abs(rmwaxis + .5 )).argmin() + len( startdist.values)
    ind_plus5 = (np.abs(rmwaxis - .5 )).argmin() + len( startdist.values)

    # find cloud heights for current eye
    H, xaxis_value = cloud_height.find_cloud_heights( crl_name, -30, ind0, ind1, xaxis ='in-situ-dist', crl_path=crl_path, new_heights=True)

    # break the distribution down by rmw distance!
    lowrmw = []
    highrmw = []
    # cycle through every case
    for i in range( len( H)):

        # high rmw case: add cloud heights
        if rmwaxis[ i] > .5 or rmwaxis[ i] < -.5:
            highrmw.append( H[ i])
        # low rmw case
        else:
            lowrmw.append( H[ i])

    # make test plots! (code below is based off scripts found in plotting/plot_rmw_chosen_eyewalls.py)
    padding = 100
    plt.figure(figsize=(18, 19))
    gridspec.GridSpec(4,6)

    title = ( "CRL Data, TC " + metadata['tc_name'] + ", "
            + metadata['dates'][ dataset] + " Eye Pass " + metadata['eye_pass'][ dataset] )

    # add histograms!
    plt.subplot2grid((4,6), (0,0), colspan=2)
    sns.histplot( y= lowrmw, kde=True, binwidth=.175, stat='probability', edgecolor='k', linewidth=2, color = 'g')
    plt.xlabel( 'Cloud Height Probability (Sums to 1)', fontsize=14)
    plt.ylabel( 'Height from Surface (Km)')
    plt.title( "Cloud Height Distribution for Inner Eye Clouds")
    plt.ylim( [-.2, 3.75])
    plt.xlim( [0, .5])
    plt.grid(True)

    axis0 = plt.gca()
    axis0.text( .4, 4.5, title, fontsize=20)


    plt.subplot2grid((4,6), (0,3), colspan=2)
    sns.histplot( y= highrmw, kde=True, binwidth=.175, stat='probability', edgecolor='k', linewidth=2, color = 'g')
    plt.xlabel( 'Cloud Height Probability', fontsize=14)
    plt.ylabel( 'Height from Surface (Km)')
    plt.title( "Cloud Height Distribution for Outer Eye Clouds")
    plt.ylim( [-.2, 3.75])
    plt.xlim( [0, .5])
    plt.grid(True)

    # plot power with original distance axis
    # plt.subplot(311)
    plt.subplot2grid((4,6), (1,0), colspan=6)
    make_plots_new_heights.plot_power_ch1(crl_path, crl_name, 'in-situ-dist')
    plt.xlim( [ - padding, padding])
    plt.ylim( [ 0, crl_data.H_max + .1])
    # add cloud height overlay
    plt.plot( xaxis_value, H, c='r', linewidth=2)

    # add a line way out of plot to add a nice red legend entry!
    plt.axvline(x=-100000, c='r', linewidth=4, label='Cloud Heights')
    plt.legend(loc='upper left', facecolor='w', framealpha=1)

    # print("power plot created")


    # plot power with new rmw axis
    # plt.subplot(312)
    plt.subplot2grid((4,6), (2,0), colspan=6)
    make_plots_new_heights.plot_T(crl_path, crl_name, 'in-situ-dist')
    plt.xlim( [ - padding, padding])
    plt.ylim( [ 0, crl_data.H_max + .1])

    # add lines denoting eye edges:
    plt.axvline( x=lim0, c='r', linewidth=4)
    plt.axvline( x=lim1, c='r', linewidth=4, label='eyewall limits')
    # add line showing newly determined tc center!
    plt.axvline( x=center_dist, c='b', linewidth=4, label='TC Center')
    # add lines showing rmw = +-.5!
    plt.axvline( x=dist[ ind_minus5], c='g', linewidth=4)
    plt.axvline( x=dist[ ind_plus5], c='g', linewidth=4, label='rmw=.5')
    plt.legend(loc='upper left', facecolor='w', framealpha=1)

    # print("temp plot created")


    # plot in situ data
    simple_flight_level_plot.plot( crl_path, crl_name, metadata, dataset, 'dist', togrid=True)
    plt.xlim( [ - padding, padding])
    plt.locator_params(axis='x', nbins=11)
    helper_fns.add_blank_colorbar()
    plt.xlabel( "Distance from TC Center (km)")

    # print("in situ plot created")
    print('Plots Created')

    if calc_stats:
        calc_stats( H, rmwaxis, xaxis_value)
        # print( 'statistics calculated')

    return H



# a helper function used to calculate statistics for cloud height distributions split
# H represents the array holding all cloud heights, rmwaxis is an array from -1 to 1
# representing the predetermined / not actually true radius of maximum wind data, and
# xaxis is the corresponding distance axis
def calc_stats( H, rmwaxis=np.nan, xaxis=np.nan, c=.95):
    print( "Number of data points:  " + str( len( H)))
    print( "Height value range:     " + str( H.min()) + " km to " + str( H.max()) + " km")

    # get rid of heights below 50m
    only_cloud_heights = H[ np.where( H > .05)[0] ]

    # print( "Height value mean (no clear air):      " + str( np.mean( only_cloud_heights)) + " km")
    # print( "Height value median:    " + str( np.median( H)) + " km")
    # print( "Standard deviation:     " + str( np.std( H)) )

    m, low, high = mean_confidence_interval( only_cloud_heights, confidence=c)
    print( "Mean confidence outputs:")
    print( "mean = " + str( m))
    print( "confidence = " + str( c))
    print( "low bound = " + str( low))
    print( "high bound = " + str( high))

    print( "Single draw confidence outputs:")
    m, low, high = single_draw_confidence( only_cloud_heights, confidence=c)
    print( "low bound = " + str( low))
    print( "high bound = " + str( high))

    # these do the same thing as the function above
    # confidence2( only_cloud_heights, confidence=c)
    # confidence3( only_cloud_heights, confidence=c)

    # check if rmw and regular axes are filled
    # true case
    if ~np.isnan( rmwaxis).all() and ~np.isnan( xaxis).all():

        # determine inner vs outer cloud heights as is done above
        lowrmw = []
        highrmw = []
        for i in range( len( H)):
            # high rmw case: add cloud heights
            if rmwaxis[ i] > .5 or rmwaxis[ i] < -.5:
                highrmw.append( H[ i])
            # low rmw case
            else:
                lowrmw.append( H[ i])

        # find mean heights for each region
        innermean = np.mean( lowrmw)
        outermean = np.mean( highrmw)

        # print out mean cloud heights for inner and outer rmw regions
        print("Inner RMW mean H:    " + str( innermean))
        print("Outer RMW mean H:    " + str( outermean))




# confidence interval function taken from stack overflow lol
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

# from the same stack overflow thread as above!
def confidence2(data, confidence=0.95):
    a = 1.0 * np.array(data)
    out = scipy.stats.t.interval( confidence, len(a)-1, loc=np.mean(a), scale=scipy.stats.sem(a))
    low_bound, high_bound = out[ 0], out[ 1]

    return np.mean( a), low_bound, high_bound

# another method for finding confidence intervals, this time taken from a different stack overflow post
# matches the two functions above
def confidence3( data, confidence=.95):
    a = 1.0 * np.array(data)
    N = len( a)
    mean, sigma = np.mean(a), np.std(a)
    out = scipy.stats.norm.interval( confidence, loc=mean, scale=sigma/math.sqrt(N))
    low_bound, high_bound = out[ 0], out[ 1]
    return mean, low_bound, high_bound

# this function finds the confidence intervals for a single draw (one random choice)
# rather than a mean draw from lots of selections. The confidence interval is way worse
# in this case :( but idk if this is the statistic that I want to use
def single_draw_confidence( data, confidence=.95):
    a = 1.0 * np.array(data)
    N = len( a)
    mean, sigma = np.mean(a), np.std(a)
    out = scipy.stats.norm.interval( confidence, loc=mean, scale=sigma)
    low_bound, high_bound = out[ 0], out[ 1]
    return mean, low_bound, high_bound




    # stull p 42 was used for this definition
    # think about nyquist frequency in relation to pulses of tc convection, stull p 306
    '''
    # the following code leads to the same output as the line above! So just calculate
    # things using numpy for speed + consistency
    # calculate and print cloud height standard deviation
    meanH = np.mean( H)
    N = len( H)
    sum = 0
    for H_i in H:
        sum += (H_i - meanH)**2
    standard_dev = ( sum / N)**.5
    print( "Calculated standard dev: " + str( standard_dev))
    '''
