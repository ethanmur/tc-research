# import...
import numpy as np
import os
import sys
import pandas as pd
from geopy import distance
import metpy.calc as mpcalc
import xarray as xr
import math
import matplotlib.pyplot as plt
os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023
sys.path.append("/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops

# Goal: create a generalized function that finds cloud heights and compiles info for all composite groups (intensity, shear, etc).
# 		do all the sorting later: the information should all be in this dataframe!
# 		can also save locally! easier access later
# 		code is partially based on "code/eye cloud paper figures/Figure 4 cloud dists vs shear quadrant.ipynb"
# Inputs: None! Just run this script to create the nice, organized data file
# 		 No need for a tc='all' input: just do the sorting after making this dataframe!
# Return: a dataframe with height / distribution information for every pass! 
# Notes * -> a list, category = WH, intensifying, DSR depending on input flags! Structure:
# flight | pass | UTC time * | x dists * | y dists * | cloud heights * | vertical bins | normalized dists | defined eyewalls or not | intensity | intensification | tc category | shear strength | shear dir
# case 0
# case 1
# ... 
def find_dists():
    # use a helper fn to get year and file names, along with intensity, etc. metadata
    tc='all'
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc=tc)
    metadata = eyewall_metadata.all_metadata()
    # the dataframe below contains simple xy distances for each day (not passes) and is locally saved. 
    # Recreate this dataset using "code/tests/2023-06-07 find cartesian and shear xy distances.ipynb"
    df_dists = pd.read_pickle("/Users/etmu9498/research/data/simple_distances.pkl")
    # use the helper function below to go from single dates to individual eye passes. Add to this dataframe for the rest of the code
    df_dists = separate_passes( df_dists)
    # add intensity, shear, category, etc info during this step!
    df_final = add_info( df_dists)
    return df_final

# after creating a dataframe with passes, cloud heights, and distances in separate columns, 
# add more info like shear and intensity here. Similar methods as before!
def add_info( df_dists):
    # the same code as above, but only do this for TC eye passes!
    tc = 'all'
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc=tc)
    metadata = eyewall_metadata.all_metadata()
    intensity_l, intensification_l, category_l, shearmag_l, sheardir_l = [], [], [], [], []

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):        
            # get eye time limits for each date
            date = crlname[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if crlname[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                    intensity = metadata[ yearval]['intensity'][ '0812am']
                    intensification = metadata[ yearval]['intensification'][ '0812am']
                    category = metadata[ yearval]['category'][ '0812am']
                    shearmag = metadata[ yearval]['shear_mag'][ '0812am']
                    sheardir = metadata[ yearval]['shear_dir'][ '0812am']
                elif crlname[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
                    intensity = metadata[ yearval]['intensity'][ '0812pm']
                    intensification = metadata[ yearval]['intensification'][ '0812pm']
                    category = metadata[ yearval]['category'][ '0812pm']
                    shearmag = metadata[ yearval]['shear_mag'][ '0812pm']
                    sheardir = metadata[ yearval]['shear_dir'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
                intensity = metadata[ yearval]['intensity'][ date]
                intensification = metadata[ yearval]['intensification'][ date]
                category = metadata[ yearval]['category'][ date]
                shearmag = metadata[ yearval]['shear_mag'][ date]
                sheardir = metadata[ yearval]['shear_dir'][ date]
       	    else:
                eyewall_limits = [ ()]        
            # do this for each of the eyewall limit pairs! Can have multiple eyes per crl dataset
            for eyei, eyeval in enumerate( eyewall_limits):
                if ~ np.isnan(intensity):
                    intensity_l.append( intensity)
                    intensification_l.append( intensification)
                    category_l.append( category)
                    shearmag_l.append( shearmag)
                    sheardir_l.append( sheardir)
    # add the new categories to the dataframe!
    df_dists['intensity'] = intensity_l
    df_dists['intensification'] = intensification_l
    df_dists['category'] = category_l
    df_dists['shearmag'] = shearmag_l
    df_dists['sheardir'] = sheardir_l
    return df_dists

# Helper fn used to take xy distances from individual days and split them into specific passes!
# inputs: df_dists, a pandas dataframe with date, time, and xy distance info for single days.
# return: df_new, a pandas dataframe with original and shear corrected xy distances for each individual eye pass.
# single layer cloud data are used for this analysis: easier to match single heights to single xy distance
def separate_passes( df_dists):
    # the same code as above, but only do this for TC eye passes!
    tc = 'all'
    multi_layers = False
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc=tc)
    metadata = eyewall_metadata.all_metadata()
    crl_root_path = "/Users/etmu9498/research/data/crl-all-data-processed/"
    # add trimmed time and distance values here
    datetrim, passtrim, timetrim, xtrim, ytrim, xsheartrim, ysheartrim, cloudheights, definedlist = [], [], [], [], [], [], [], [], []
    alltimes = df_dists['times'].values
    allxdists = df_dists['xdists'].values
    allydists = df_dists['ydists'].values
    shearxdists, shearydists = find_shear_dists(df_dists) # use helper fn to find shear shifted xy distances

    # keep track of which tc we're on with this counter
    flightcounter = 0
    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):    
            # get eye time limits for each date
            date = crlname[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if crlname[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                    defined_eyewalls = metadata[ yearval]['defined_eyewall'][ '0812am']
                elif crlname[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
                    defined_eyewalls = metadata[ yearval]['defined_eyewall'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
                defined_eyewalls = metadata[ yearval]['defined_eyewall'][ date]
            else:
                eyewall_limits = [ ()]
            # do this for each of the eyewall limit pairs! Can have multiple eyes per crl dataset
            for eyei, eyeval in enumerate( eyewall_limits):
                if len( eyeval) > 0:
                    # find the corresponding indices to the time limits
                    ind0 = np.argmin( np.abs( alltimes[flightcounter] - eyeval[0] ))
                    ind1 = np.argmin( np.abs( alltimes[flightcounter] - eyeval[1] ))
                    # clip relevant fields down to the eyewall limits
                    # load crl data
                    os.chdir( crl_root_path + yearval)
                    crl_data = xr.open_dataset( crlname)
                    H = crl_data.height
                    power = crl_data.P_ch1[ ind0 : ind1, :]
                    axis = crl_data.time[ ind0 : ind1]
                    p3_height = crl_data.p3_height[ ind0 : ind1]
                    # find cloud top heights for values within the specified eye distance range
                    if yearval == '2021':
                        cutoff = -30
                    elif yearval == '2022':
                        cutoff = -40
                    if multi_layers:
                        heights_lists, time, count = find_cloud_tops.find_multi_cloud_heights( H, power, axis, p3_height, cutoff_power = cutoff)
                        heights = np.array( [item for sublist in heights_lists for item in sublist])
                    else:
                        # regular case: find the top cloud height layer
                        heights, time = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = cutoff)
    				# add xy and height data to the proper lists	                
                    if ind0 != ind1:
                        datetrim.append( crlname[3:13])
                        passtrim.append( eyei)
                        xtrim.append( allxdists[flightcounter][ind0:ind1])
                        ytrim.append( allydists[flightcounter][ind0:ind1])
                        xsheartrim.append( shearxdists[flightcounter][ind0:ind1])
                        ysheartrim.append( shearydists[flightcounter][ind0:ind1])
                        cloudheights.append( heights)
                        timetrim.append( alltimes[flightcounter][ind0:ind1])
                        definedlist.append( defined_eyewalls)        
                    else:
                        print("Error: matching indices :(")
            # skip over two dates that create problems for some reason
            if crlname == 'P3_20220831H1_processed.nc' or crlname == 'P3_20220830H1_processed.nc':
                pass # print('ERROR CASE!!')
            else:
                flightcounter += 1
    # convert lists to a nice dataframe!
    df_dists_trim = pd.DataFrame()
    df_dists_trim['flight'] = datetrim
    df_dists_trim['pass'] = passtrim
    df_dists_trim['times'] = timetrim
    df_dists_trim['xdists'] = xtrim
    df_dists_trim['ydists'] = ytrim
    df_dists_trim['xdistsshear'] = xsheartrim
    df_dists_trim['ydistsshear'] = ysheartrim
    df_dists_trim['cloudheights'] = cloudheights
    df_dists_trim['Defined Eyewalls'] = definedlist
    return df_dists_trim

# helper fns (x and y are arrays)
def cart2pol(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return(r, theta)
def pol2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return(x, y)

# go from cartesian xy to shear shifted xy distances here
def find_shear_dists( df_dists):
    yearlist, crlfilelist = helper_fns_winter2023.get_crl_datasets( tc='all')
    metadata = eyewall_metadata.all_metadata()
    allxdists = df_dists['xdists'].values
    allydists = df_dists['ydists'].values

    sheardirlist = []
    # pulling shear for each case- make sure not to get shear from 08/30 and 08/31 which have no data!
    # keep track of which tc we're on with this counter
    flightcounter = 0
    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        for filei, crlname in enumerate( crlfilelist[ yeari]):        
            # get eye time limits for each date
            date = crlname[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if crlname[11:13] == "H1":
                    sheardir = metadata[ yearval]['shear_dir'][ '0812am']
                elif crlname[11:13] == "H2":
                    sheardir = metadata[ yearval]['shear_dir'][ '0812pm']
            elif date in metadata[ yearval]['shear_dir'].keys():
                sheardir = metadata[ yearval]['shear_dir'][ date]
            else:
                sheardir = [ ()]	        
            # deal with problem cases
            if crlname == 'P3_20220831H1_processed.nc' or crlname == 'P3_20220830H1_processed.nc':
                pass # print('ERROR CASE!!')
            else:
                sheardirlist.append( sheardir)
                flightcounter += 1
    # save new x and y dists below
    shearxdists, shearydists = [], []
    # do this for every set of x and y value
    for i in range(len(sheardirlist)):
        r, theta = cart2pol( np.array(allxdists[i]), np.array(allydists[i]))	    
        theta2 = theta + np.radians( sheardirlist[i])	    
        xtemp, ytemp = pol2cart(r, theta2)
        shearxdists.append(xtemp)
        shearydists.append(ytemp)
    return shearxdists, shearydists







################################################################
# Goal: create one function to plot all intensity, intensification, and shear distributions consistently!
#       It should be much more flexible and modular.
# Inputs: tc- can be 'all' to use all cases, '2021' or '2022' for one year, 'defined' to just use defined eyewalls, 
#         or a list of indices for specific passes.
#         bintype- whether to make a 'category', 'intensity', 'intensification', 'shearmag', or 'sheardir' composite plot.
#         newlevels- defaults to false, and special levels are used depending on bintype. Otherwise, give a list of 
#         additional values / levels to sort by.
# Return: None. Side effect: create composite distributions.
def plot_composites( tc='all', df_passes=False, bintype='category', newlevels=False, newlabels=False, plottype='averaged', xaxistype='probability-norm', binwidth=25, smoothwidth=25, figtitle=False):
    # first, load the whole set of eye passes here
    # if already created / passed to this function, continue. otherwise, make the dataset.
    if not isinstance(df_passes, pd.DataFrame):
        df_passes = find_dists()
    print("All eye passes loaded.")

    # next, depending on tc type, remove certain eye passes from analysis
    df_passes = sorting_helper(tc, df_passes)

    # create different labels, etc depending on the plot type!
    if bintype == 'category':
        default_levels = ['td', 'ts', 'wh', 'sh']
        default_labels = ['TD', 'TS', 'WH', 'SH']
        title = "Category"
    elif bintype == 'intensity':
        default_levels = [50, 100]
        default_labels = ['Weak', 'Moderate', 'Strong']
        title = "Intensity"
    elif bintype == 'intensification':
        default_levels = [-7.5, 7.5]
        default_labels = ['Weakening', 'Steady State', 'Strengthening']
        title = 'Intensification'
    elif bintype == 'shearmag':
        default_levels = [ 7.5, 15]
        default_labels = ['Weak', 'Moderate', 'Strong']
        title = "Shear Magnitude"
    elif bintype == 'sheardir':
        default_levels = ['DR', 'DL', 'UL', 'UR']
        default_labels = default_levels
        title = "Shear Quadrant"
        
    # choose to use the default or manually entered levels / composite bins here
    if newlevels:
        levels = newlevels
    else:
        levels = default_levels
    if newlabels:
        labels = newlabels
    else:
        labels = default_labels
    
    # save / concat cloud heights here! this list of lists is the length of the levels being used (1 list per category).
    # cloud_composites will hold separate lists for each pass, 
    # while cloud_composites_flat will flatten all heights into one big list.
    cloud_composites, cloud_composites_flat = [], []

    # normal case: just add 4 labels
    if bintype == 'category' or bintype == 'sheardir':
        for i in range(len(levels)):
            cloud_composites.append([])
            cloud_composites_flat.append([])
    # bin value cases: add an extra list for above / below bounds!
    else:
        for i in range(len(levels) + 1):
            cloud_composites.append([])
            cloud_composites_flat.append([])
    # Cycle through eye passes here
    for passi, passval in enumerate( df_passes['pass']):
        # do the sorting here! the special cases are taken care of first, then the default cases
        # (works for intensity, intensification, shearmag) come next
        
        # category case- find the index for the matching category and save cloud heights there
        if bintype == 'category':
            cat = df_passes['category'][passi] # category for this pass
            catlisti = levels.index(cat) # the index in the default or given levels list (where to append cloud heights)
            cloud_composites[catlisti].append(df_passes['cloudheights'][passi].tolist()) # add cloud heights!
            cloud_composites_flat[catlisti] += df_passes['cloudheights'][passi].tolist()

        # shear direction case: relatively complicated. Need to search through each cloud height, not entire passes
        elif bintype == 'sheardir':
            # save temp shear dir values here for each pass
            tempshears = [ [], [], [], []]
            # go through each point in each eye pass
            for cloudi in range(len( df_passes['cloudheights'][passi])):
                # current distances
                xi, yi = df_passes['xdistsshear'][passi][cloudi], df_passes['ydistsshear'][passi][cloudi]
                # DL case
                if xi < 0.0 and yi > 0.0:
                    datacase = 'DL'
                    quadrantlisti = levels.index(datacase) # get index of DL in list
                    tempshears[quadrantlisti] += [df_passes['cloudheights'][passi][cloudi]] # .tolist()
                # DR case
                elif xi > 0.0 and yi > 0.0:
                    datacase = 'DR'
                    quadrantlisti = levels.index(datacase)
                    tempshears[quadrantlisti] += [df_passes['cloudheights'][passi][cloudi]]
                # UR case
                elif xi > 0.0 and yi < 0.0:
                    datacase = 'UR'
                    quadrantlisti = levels.index(datacase)
                    tempshears[quadrantlisti] += [df_passes['cloudheights'][passi][cloudi]]
                # UL case
                elif xi < 0.0 and yi < 0.0:
                    datacase = 'UL'
                    quadrantlisti = levels.index(datacase)
                    tempshears[quadrantlisti] += [df_passes['cloudheights'][passi][cloudi]]
            # after looping through all cloud heights, add shear quads for this pass to cloud composite lists
            for tempsheari in range(len(tempshears)):
                cloud_composites_flat[tempsheari] += tempshears[tempsheari]
                if len(tempshears[tempsheari]) > 0:
                    cloud_composites[tempsheari].append(tempshears[tempsheari])

        # intensity / shear mag / intensification case:
        else:
            if bintype == 'intensity':
                intensity = df_passes['intensity'][passi] 
            elif bintype == 'intensification':
                intensity = df_passes['intensification'][passi]
            elif bintype == 'shearmag':
                intensity = df_passes['shearmag'][passi]
            # figure out which level the intensity level falls between.
            # if 3 numbers are present, create 4 levels (lower / upper bounds, 2 middle levels)
            # ex: | 40 kt | 80 kt | 120 kt |
            # lowest level case:
            if intensity < levels[0]:
                intensitylisti = 0
            # highest level case:
            elif intensity >= levels[-1]:
                intensitylisti = len(levels)
            # middle cases: iterate through the indices to find a matching case!
            else:
                for intensityi in range(len(levels) - 1):
                    currentval, nextval = levels[intensityi], levels[intensityi + 1]
                    if intensity >= currentval and intensity < nextval:
                        intensitylisti = 1 + intensityi
                        break
            cloud_composites[intensitylisti].append(df_passes['cloudheights'][passi].tolist()) # add cloud heights!
            cloud_composites_flat[intensitylisti] += df_passes['cloudheights'][passi].tolist()


    # after cycling through each pass, print out helpful statistics and choose the type of plot to create!
    simple_stats(cloud_composites, cloud_composites_flat, labels, bintype)
    if plottype == 'individual':
        if figtitle:
            plot_individual_dists(cloud_composites, figtitle, labels, binwidth, smoothwidth, xaxistype)
        else:
            plot_individual_dists(cloud_composites, title, labels, binwidth, smoothwidth, xaxistype)

    elif plottype == 'averaged':
        if figtitle:
            plot_averages(cloud_composites_flat, figtitle, labels, binwidth, smoothwidth, xaxistype)
        else:
            plot_averages(cloud_composites_flat, title, labels, binwidth, smoothwidth, xaxistype)
        
    elif plottype == 'all':
        if figtitle:
            plot_individual_dists(cloud_composites, figtitle, labels, binwidth, smoothwidth, xaxistype)
            plot_averages(cloud_composites_flat, figtitle, labels, binwidth, smoothwidth, xaxistype)
        else:
            plot_individual_dists(cloud_composites, title, labels, binwidth, smoothwidth, xaxistype)
            plot_averages(cloud_composites_flat, title, labels, binwidth, smoothwidth, xaxistype)

    return df_passes

# a helper function used to quickly sort / remove certain cases depending on the TC flag
def sorting_helper(tc, df_passes):
    # 'all' case: keep all passes
    if tc == 'all':
        pass
    # '2021' or '2022' year case: only keep this one year
    elif tc == '2021' or tc == '2022':
        df_passes = df_passes.drop(index=np.where( ~ (df_passes['flight'].str[0:4] == tc ))[0])
        df_passes = df_passes.reset_index()
    # defined eyewall case: only keep cases with clear eyewall limits
    elif tc == 'defined':
        df_passes = df_passes.drop(index=np.where( ~ df_passes['Defined Eyewalls'])[0])
        df_passes = df_passes.reset_index()
    # list index case: drop cases at specified indices!
    elif isinstance(tc, list):
        df_passes = df_passes.drop(index=np.array(tc))
        try:
            df_passes = df_passes.reset_index()
        except:
            print("value error :/")
            print(df_passes)
    else:
        print("Please enter a valid TC type!")
        return
    return df_passes


def simple_stats(cloud_composites, cloud_composites_flat, labels, bintype):
    totalsum = 0
    for category in cloud_composites_flat:
        totalsum += len(category)
    print("Total number of lidar data points: " + str(totalsum))
    
    for datai, dataval in enumerate( cloud_composites_flat):
        test_ht = np.array( dataval)
        test_ht = test_ht[ np.where( test_ht > 50)[0] ]
        print( "\n" + labels[ datai])
        # don't print cases for shear dir data- doesn't really make sense
        if bintype != 'sheardir':
            print( "cases: " + str( len( cloud_composites[datai]) ))
        print("orig number of data points: " + str( len( dataval)))
        print("non zero number of data points: " + str( len( test_ht)))
        print("orig mean ht = " + str( np.round( np.nanmean( dataval), 3)))
        print("non zero mean ht = " + str( np.round( np.nanmean( test_ht), 3)))
        if len( dataval) > 0:
            print("clear air % = " + str( 100 *( 1 - np.round( len( test_ht ) / len( dataval), 3)) ))
        else:
            print("No cases for this intensity and year combination.")


def plot_individual_dists(cloud_composites, title, labels, binwidth, smoothwidth, xaxistype, lw=2):
    colors = ['b', 'k', 'y', 'g', 'c', 'Grey']

    # do this for each of the composite categories:
    for categoryi in range(len( cloud_composites)):

        # make figure before loop
        fig, a0 = plt.subplots(1, 1, figsize=( 5, 5) )
        fs = 16
        # loop through each case

        print("Category " + labels[categoryi])
        for heighti, cloudheight in enumerate(cloud_composites[categoryi]):
            # remove 0 km heights from figure and define a cloud bin axis
            origheight = np.array(cloudheight)
            cloudheight = origheight[ np.where( origheight > 50)[0] ]
            height_bin=np.arange(0, 4510, binwidth)
            # use a helper function to bin, smooth, and normalize data
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)

            probs_sum = 0
            for height in height_probs:
                probs_sum += height
            # print("Case " + str(heighti) + ": sums to " + str(np.round(probs_sum, 2)) + "%, cloudy air = " + str( (len(cloudheight) / len(origheight)) * 100) + " %")
            
            if xaxistype == 'probability':
                # add the probability plot to the total figure for this category
                a0.plot( height_probs, heights, color=colors[categoryi], linewidth=lw)
            elif xaxistype == 'probability-norm':
                # add the probability plot to the total figure for this category
                a0.plot( height_probs_norm, heights_norm, color=colors[categoryi], linewidth=lw)    
       
        # save the total figure!
        a0.set_ylabel( 'Height from Surface (Km)', fontsize=fs)
        if xaxistype == 'probability' or xaxistype=='probability-norm':
            a0.set_xlabel('Cloud Height Probability (%)', fontsize=fs)
        a0.set_ylim( [-200, 4600])
        a0.grid(False)
        a0.set_title( "Cloud Heights vs TC " + title + ", " + labels[categoryi] + " Case", fontsize=fs)    
        a0.tick_params(axis='both', which='major', labelsize=fs)

def plot_averages(cloud_composites_flat, title, labels, binwidth, smoothwidth, xaxistype, lw=2):
    colors = ['b', 'k', 'y', 'g', 'c', 'Grey']
    # make figure before loop
    fig, a0 = plt.subplots(1, 1, figsize=( 7, 7) )
    fs = 16
    # loop through each case
    for heighti, cloudheight in enumerate(cloud_composites_flat):
        # remove 0 km heights from figure and define a cloud bin axis
        origheight = np.array(cloudheight)
        cloudheight = origheight[ np.where( origheight > 50)[0] ]
        height_bin=np.arange(0, 4510, binwidth)
        # use a helper function to bin, smooth, and normalize data
        

        if xaxistype == 'probability':
            # height_probs, heights = binning_orig_method( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            # add the probability plot to the total figure! will add nice touches and save later
            a0.plot( height_probs, heights, color=colors[heighti], label=labels[heighti], linewidth=lw)
            if heighti == 0:
                a0.set_xlabel('Cloud Height Probability (%)', fontsize=fs)
        elif xaxistype == 'probability-norm':
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            # add the probability plot to the total figure! will add nice touches and save later
            a0.plot( height_probs_norm, heights_norm, color=colors[heighti], label=labels[heighti], linewidth=lw)
            if heighti == 0:
                a0.set_xlabel('Cloud Height Probability (%)', fontsize=fs)
    
    a0.set_ylabel( 'Height from Surface (Km)', fontsize=fs)
    a0.set_ylim( [-200, 4600])
    a0.grid(False)
    a0.set_title( "Cloud Heights vs TC " + title, fontsize=fs)    
    a0.tick_params(axis='both', which='major', labelsize=fs)
    plt.legend(fontsize=fs, loc="upper right")


def binning_orig_method( origheight, cloudheight, height_bin, binwidth, smoothwidth):
    # calculate mean height variation within the bin range
    mean_height=[]
    mean_height_count=[]
    mean_height_prob = []
    # do this for every height level determined by the manually inputed bin width
    for newi in height_bin:
        # find the points that fall within this height bin for this step. res = indices of valid heights
        res=np.where(np.logical_and( cloudheight >= newi - binwidth / 2., cloudheight <= newi + binwidth / 2. ))
        # add mean heights of points in this bin, and the count, to respective lists
        mean_height.append( np.mean( cloudheight[ res]))
        mean_height_count.append( len( res[0]))
        # use origheight, not cloudheight, to include heights below 50m in normalization.
        # this line accounts for clear air fraction when scaling the curves!
        if len( origheight) > 0:
            mean_height_prob.append( len( res[0] ) / len( origheight))

        # this should only happen if there's no data at this height bin... unlikely,
        # unless there are no cases at this cat
        else:
            # just give this bin a really high height?
            # mean_height_prob.append( len( res[0] ) / .00001 )
            # or, just make it 1 (all points here)?
            mean_height_prob.append(1)

    # smooth data before plotting to eliminate noise
    box_pts = smoothwidth
    box = np.ones(box_pts)/box_pts
    prob_smooth = np.convolve( mean_height_prob, box, mode='same')

    # remove 0% probability lines (mostly at high alts) from histogram plots!
    prob_smooth_trim = prob_smooth[ np.where( prob_smooth > 0.00001)[0]]
    height_trim = height_bin[ np.where( prob_smooth > 0.00001)[0]]

    return prob_smooth_trim * 100, height_trim


# mostly the same as the code above, but it replaces the confusing convolve function with a simple boxcar.
# results are the same luckily, no matter the smoothing function!
def binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth):
    # calculate mean height variation within the bin range
    mean_height=[] # good for testing bin sorting, but not returned
    mean_height_count=[] # not sure why i'm doing this lol
    mean_height_prob = []
    mean_height_prob_normalized = []

    # do this for every height level determined by the manually inputed bin width
    for newi in height_bin:
        # find the points that fall within this height bin for this step. res = indices of valid heights
        res=np.where(np.logical_and( cloudheight >= newi - binwidth / 2., cloudheight <= newi + binwidth / 2. ))
        # add mean heights of points in this bin, and the count, to respective lists
        mean_height.append( np.mean( cloudheight[ res]))
        mean_height_count.append( len( res[0]))

        # use origheight, not cloudheight, to include heights below 50m in normalization.
        # this line accounts for clear air fraction when scaling the curves!
        if len( origheight) > 0:
            mean_height_prob.append( len( res[0] ) / len( cloudheight))
            mean_height_prob_normalized.append( len( res[0] ) / len( origheight))

        # this should only happen if there's no data at this height bin... unlikely,
        # unless there are no cases at this cat
        else:
            # just give this bin a really high height?
            mean_height_prob.append( len( res[0] ) / .00001 )
            mean_height_prob_normalized.append( len( res[0] ) / .00001 )
            # or, just make it 1 (all points here)?
            # mean_height_prob.append(1)
            # mean_height_prob_normalized.append(1)

    # smooth data before plotting to eliminate noise
    mean_height_prob = np.array(pd.Series( np.array(mean_height_prob)).rolling(window = smoothwidth, min_periods=1, center=True).mean())
    mean_height_prob_normalized = np.array(pd.Series( np.array(mean_height_prob_normalized)).rolling(window = smoothwidth, min_periods=1, center=True).mean())

    # remove 0% probability lines (mostly at high alts) from histogram plots!
    height_trim = height_bin[ np.where( mean_height_prob > 0.00001)[0]]
    height_norm = height_bin[ np.where( mean_height_prob_normalized > 0.00001)[0]]
    mean_height_prob = mean_height_prob[ np.where( mean_height_prob > 0.00001)[0]]
    mean_height_prob_normalized = mean_height_prob_normalized[ np.where( mean_height_prob_normalized > 0.00001)[0]]

    return mean_height_prob * 100, mean_height_prob_normalized * 100, height_trim, height_norm
















# This script focuses on creating monte carlo distributions for a specific case (TDs, strengthening systems, etc).
# All the sorting has been done in a previous script! So, skip those steps. Also, create a new monte carlo figure with all dists plotted atop one another
def plot_monte_carlo( keepinds, removeinds, df_passes=False, newlabels=False, xaxistype='probability-norm', binwidth=25, smoothwidth=25, figtitle=False, color='b', title=''):
    # first, load the whole set of eye passes here if not passed as an arguement
    if not isinstance(df_passes, pd.DataFrame):
        df_passes = find_dists()          
    # plot setup  
    if newlabels:
        labels = newlabels
    else:
        labels = "monte carlo test"
    
    # save the distribution with all all cases (the original distribution plot) and a list of monte carlo passes in the following lists:
    # make sure the distributions are flat for easier plotting.
    orig_dist, monte_dists = [], []

    # create the original distirbution here
    orig_passes = df_passes.drop(index=np.array(removeinds))
    orig_passes = orig_passes.reset_index()
    for passi, passval in enumerate( orig_passes['pass']):
            orig_dist += orig_passes['cloudheights'][passi].tolist()

 
    # Cycle through eye passes here
    for keepind in keepinds:
        # remove all removeinds, along with the current keepind, for this current monte carlo test.
        currentremove = np.append(removeinds, keepind)  
        df_passes_current = df_passes.drop(index=np.array(currentremove))
        df_passes_current = df_passes_current.reset_index()

        # find cloud heights for remaining inds
        monte_dists_i = [] # store each pass for this iteration here. append after.
        for passi, passval in enumerate( df_passes_current['pass']):
            monte_dists_i += df_passes_current['cloudheights'][passi].tolist()
        monte_dists.append( monte_dists_i)


    print("Number of data points in original dist: " + str( len(orig_dist)))
    print("Number of Monte Carlo Tests: " + str(len(monte_dists)))
    # print(orig_dist)

    #########
    # make the distributions for each case (original and monte carlos) and plot them!!
    #########

    # make figure before loop
    fig, a0 = plt.subplots(1, 1, figsize=( 7, 7) )
    fs = 16
    lw = 2
    a = .5
    # loop through each monte carlo case
    for heighti, cloudheight in enumerate(monte_dists):

        # print("Height count for case " + str(heighti) + ": " + str(len( cloudheight)))
        
        # remove 0 km heights from figure and define a cloud bin axis
        origheight = np.array(cloudheight)
        cloudheight = origheight[ np.where( origheight > 50)[0] ]
        height_bin=np.arange(0, 4510, binwidth)
        # use a helper function to bin, smooth, and normalize data
        
        if xaxistype == 'probability':
            # height_probs, heights = binning_orig_method( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
        
            # add the probability plot to the total figure! will add nice touches and save later
            if heighti == 0:
                a0.plot( height_probs, heights, color=color, label='Monte Carlos', linewidth=lw, alpha=a)
            else:
                a0.plot( height_probs, heights, color=color, linewidth=lw, alpha=a)
            
        elif xaxistype == 'probability-norm':
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            # add the probability plot to the total figure! will add nice touches and save later
            if heighti == 0:
                a0.plot( height_probs_norm, heights_norm, color=color, label='Monte Carlos', linewidth=lw, alpha=a)
            else:
                a0.plot( height_probs_norm, heights_norm, color=color, linewidth=lw, alpha=a)


    # plot the original distribution last here! same steps as above
    origheight = np.array(orig_dist)
    cloudheight = origheight[ np.where( origheight > 50)[0] ]
    height_bin=np.arange(0, 4510, binwidth)    
    if xaxistype == 'probability':
        # height_probs, heights = binning_orig_method( origheight, cloudheight, height_bin, binwidth, smoothwidth)
        height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
    
        # add the probability plot to the total figure! will add nice touches and save later
        a0.plot( height_probs, heights, color=color, label='Original Dist', linewidth=lw, alpha=1)
    elif xaxistype == 'probability-norm':
        height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
        # add the probability plot to the total figure! will add nice touches and save later
        a0.plot( height_probs_norm, heights_norm, color=color, label='Original Dist', linewidth=lw, alpha=1)

    a0.set_xlabel('Cloud Height Probability (%)', fontsize=fs)
    a0.set_ylabel( 'Height from Surface (Km)', fontsize=fs)
    a0.set_ylim( [-200, 4600])
    a0.grid(False)
    a0.set_title( "Jackknife Test" + title, fontsize=fs)    
    a0.tick_params(axis='both', which='major', labelsize=fs)
    # plt.legend(fontsize=fs/2, loc="upper right", bbox_to_anchor=( 1.32, 1))
    plt.legend(fontsize=fs, loc="upper right")



# Much like the first monte carlo function, but intead of removing single passes, it removes random data points.
# nice because this test can be repeated thousands of times to give a confidence interval?!
# removecount is the number of random height points to remove. removeinds are the indices of passes to remove (to sort for TDs, strong shear, etc)
def plot_monte_carlo_random_remove( removecount, ntests, removeinds, df_passes=False, newlabels=False, xaxistype='probability-norm', binwidth=25, smoothwidth=25, figtitle=False, color='b', title=''):
    # first, load the whole set of eye passes here if not passed as an arguement
    if not isinstance(df_passes, pd.DataFrame):
        df_passes = find_dists()
    print("All eye passes loaded.")
          
    # plot setup  
    if newlabels:
        labels = newlabels
    else:
        labels = "monte carlo test"
    
    # save the distribution with all all cases (the original distribution plot) and a list of monte carlo passes in the following lists:
    # make sure the distributions are flat for easier plotting.
    orig_dist, monte_dists = [], []

    # create the original distirbution here
    orig_passes = df_passes.drop(index=np.array(removeinds))
    orig_passes = orig_passes.reset_index()
    for passi, passval in enumerate( orig_passes['pass']):
            orig_dist += orig_passes['cloudheights'][passi].tolist()

 
    # Cycle through eye passes here. 
    for keepind in range(ntests):
        # remove removecount number of random heights here. Append this list to the monte_dists list.
        randinds = np.random.randint( 0, len(orig_dist), size=removecount)
        monte_dists_i = np.delete( np.array( orig_dist), randinds)
        monte_dists.append( monte_dists_i.tolist())


    print("Number of data points in original dist: " + str( len(orig_dist)))
    print("Number of Monte Carlo Tests: " + str(len(monte_dists)))
    # print(orig_dist)

    #########
    # make the distributions for each case (original and monte carlos) and plot them!!
    #########

    # make figure before loop
    fig, a0 = plt.subplots(1, 1, figsize=( 7, 7) )
    fs = 16
    lw = 2
    a = .25
    # loop through each monte carlo case
    for heighti, cloudheight in enumerate(monte_dists):

        # print("Height count for case " + str(heighti) + ": " + str(len( cloudheight)))
        
        # remove 0 km heights from figure and define a cloud bin axis
        origheight = np.array(cloudheight)
        cloudheight = origheight[ np.where( origheight > 50)[0] ]
        height_bin=np.arange(0, 4510, binwidth)
        # use a helper function to bin, smooth, and normalize data
        
        if xaxistype == 'probability':
            # height_probs, heights = binning_orig_method( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)        
            # add the probability plot to the total figure! will add nice touches and save later
            if heighti == 0:
                a0.plot( height_probs, heights, color=color, label='Monte Carlos', linewidth=lw, alpha=a)
            else:
                a0.plot( height_probs, heights, color=color, linewidth=lw, alpha=a)
            
        elif xaxistype == 'probability-norm':
            height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
            # add the probability plot to the total figure! will add nice touches and save later
            if heighti == 0:
                a0.plot( height_probs_norm, heights_norm, color=color, label='Monte Carlos', linewidth=lw, alpha=a)
            else:
                a0.plot( height_probs_norm, heights_norm, color=color, linewidth=lw, alpha=a)


    # fancy color selection
    if color == 'k':
        c = 'c'
    else:
        c = 'k'

    # plot the original distribution last here! same steps as above
    origheight = np.array(orig_dist)
    cloudheight = origheight[ np.where( origheight > 50)[0] ]
    height_bin=np.arange(0, 4510, binwidth)    
    if xaxistype == 'probability':
        # height_probs, heights = binning_orig_method( origheight, cloudheight, height_bin, binwidth, smoothwidth)
        height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
    
        # add the probability plot to the total figure! will add nice touches and save later
        a0.plot( height_probs, heights, color=c, label='Original Dist', linewidth=lw, alpha=1)
    elif xaxistype == 'probability-norm':
        height_probs, height_probs_norm, heights, heights_norm = binning_boxcar( origheight, cloudheight, height_bin, binwidth, smoothwidth)
        # add the probability plot to the total figure! will add nice touches and save later
        a0.plot( height_probs_norm, heights_norm, color=c, label='Original Dist', linewidth=lw, alpha=1)

    a0.set_xlabel('Cloud Height Probability (%)', fontsize=fs)
    a0.set_ylabel( 'Height from Surface (Km)', fontsize=fs)
    a0.set_ylim( [-200, 4600])
    a0.grid(False)
    a0.set_title( "Monte Carlo Test" + title, fontsize=fs)    
    a0.tick_params(axis='both', which='major', labelsize=fs)
    # plt.legend(fontsize=fs/2, loc="upper right", bbox_to_anchor=( 1.32, 1))
    plt.legend(fontsize=fs, loc="upper right")