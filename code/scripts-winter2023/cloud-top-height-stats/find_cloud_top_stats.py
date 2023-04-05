import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import matplotlib.ticker as ticker


import xarray as xr
import seaborn as sns # for making prettier pdfs
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import helper_fns
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023/")
import helper_fns_winter2023
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023/cloud-top-height-stats")
import eyewall_metadata
import find_cloud_tops
import cloud_top_plotting


# make eye cloud height distributions for the four tc categories! see how they compare.
# confirm that results match earlier 2021 data, and repeat tests for 2022 data.
# much of this code is taken from "scripts/statistics/cloud_height_pdfs_all_one_figure.py/cloud_height_vs_intensity()"
def find_intensity_stats( tc='all', binwidth=25, smoothwidth=25, eye_limits='default'):
    # empty lists that will hold all the height datasets for each intensity category
    td_heights, td_cases = [], 0
    ts_heights, ts_cases = [], 0
    wh_heights, wh_cases = [], 0
    sh_heights, sh_cases = [], 0
    crl_root_path = "/Users/etmu9498/research/data/crl-all-data-processed/"


    # use a helper fn to get the relevant years and files
    yearlist, filelist = helper_fns_winter2023.get_crl_datasets( tc=tc)

    # print out the number of files to be used
    filecount = 0
    for yeari in range( len( filelist)):
        # count all the names in this year, and add to the count
        filecount += len( filelist[ yeari])
    print("Number of data files to be used in analysis: " + str( filecount))

    # load eyewall limits from helper function
    metadata = eyewall_metadata.all_metadata( eye_limits=eye_limits)

    # do this for all the datasets! years and filenames
    for yeari, yearval in enumerate( yearlist):
        print( "Year: " + yearval)
        for filei, fileval in enumerate( filelist[ yeari]):
            print( "File: " + fileval)

            # grab the limits for this case
            # simplified filename
            date = fileval[7:11]
            # check if this date exists... if not, give it some empty eyewall limits!
            # also account for fred am and pm cases!!
            if date == '0812':
                if fileval[11:13] == "H1":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812am']
                elif fileval[11:13] == "H2":
                    eyewall_limits = metadata[ yearval]['eyewall_limits'][ '0812pm']
            elif date in metadata[ yearval]['eyewall_limits'].keys():
                eyewall_limits = metadata[ yearval]['eyewall_limits'][ date]
            else:
                eyewall_limits = [ ()]

            # make sure intensity metadata are inputed for this year!
            if len( metadata[yearval]['intensity'].keys() ) > 0:
                if date == '0812':
                    if fileval[11:13] == "H1":
                        intensity = metadata[ yearval]['intensity'][ '0812am']
                    elif fileval[11:13] == "H2":
                        intensity = metadata[ yearval]['intensity'][ '0812pm']
                elif date in metadata[ yearval]['intensity'].keys():
                    intensity = metadata[yearval]['intensity'][date]
                else:
                    print( metadata[ yearval]['intensity'].keys())
                    intensity = 0
            else:
                intensity = 0

            # do this for each of the eyewall limit pairs! Can have multiple eyes per crl dataset
            for eyei, eyeval in enumerate( eyewall_limits):
                # load crl data
                os.chdir( crl_root_path + yearval)
                crl_data = xr.open_dataset( fileval)

                if len( eyeval) > 0:

                    # find the corresponding indices to the time limits
                    ind0 = np.argmin( np.abs(crl_data.time.values - eyeval[0] ))
                    ind1 = np.argmin( np.abs(crl_data.time.values - eyeval[1] ))

                    # clip relevant fields down to the eyewall limits
                    H = crl_data.height
                    power = crl_data.P_ch1[ ind0 : ind1, :]
                    axis = crl_data.time[ ind0 : ind1]
                    p3_height = crl_data.p3_height[ ind0 : ind1]

                    #print( "h: " + str( len( H)))
                    #print( "power: " + str( np.shape( power)))
                    #print( 't: ' + str( len( axis)))

                    # find cloud top heights for values within the specified eye distance range
                    if yearval == '2021':
                        cutoff = -30
                    elif yearval == '2022':
                        cutoff = -40
                    heights, time = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = cutoff)

                    # save height values from this run
                    if date == '0812':
                        if fileval[11:13] == "H1":
                            category = metadata[yearval]['category']['0812am'] # tc intensity category
                        elif fileval[11:13] == "H2":
                            category = metadata[yearval]['category']['0812pm'] # tc intensity category
                    elif date in metadata[ yearval]['eyewall_limits'].keys():
                        category = metadata[yearval]['category'][date]

                    # figure out where to put the height data depending on the tc intensity!
                    if category == 'td':
                        td_heights += np.ndarray.tolist( heights)
                        # print( 'td')
                        td_cases += 1
                    elif category == 'ts':
                        ts_heights += np.ndarray.tolist( heights)
                        ts_cases += 1
                    elif category == 'wh':
                        wh_heights += np.ndarray.tolist( heights)
                        wh_cases += 1
                    elif category == 'sh':
                        sh_heights += np.ndarray.tolist( heights)
                        sh_cases += 1

                    # optional: make a simple plot showing determined eye cloud heights!
                    #cloud_top_plotting.eye_cloud_tops( heights, time, axis, H, power, title=fileval + '-eye-' + str( eyei))
                    #os.chdir( "/Users/etmu9498/research/figures/CRL-all-data-processed/2021-eye-clouds/")
                    #plt.savefig( fileval + '-eye-' + str( eyei) + ".png", dpi = 250)


    data = [ td_heights, ts_heights, wh_heights, sh_heights]
    cases = [ td_cases, ts_cases, wh_cases, sh_cases]
    tccat = ['td', 'ts', 'wh', 'sh']
    for datai, dataval in enumerate( data):
        test_ht = np.array( dataval)
        test_ht = test_ht[ np.where( test_ht > 50)[0] ]

        print( "\n" + tccat[ datai])
        print( "cases: " + str( cases[ datai]))
        print("orig number of data points: " + str( len( dataval)))
        print("non zero number of data points: " + str( len( test_ht)))
        print("orig mean ht = " + str( np.round( np.nanmean( dataval), 3)))
        print("non zero mean ht = " + str( np.round( np.nanmean( test_ht), 3)))
        if len( dataval) > 0:
            print("clear air % = " + str( 100 *( 1 - np.round( len( test_ht ) / len( dataval), 3)) ))
        else:
            print("No cases for this intensity and year combination.")

    # make a nice figure for these tc intensities!
    cloud_top_plotting.intensity_bins( tc, td_heights, ts_heights, wh_heights, sh_heights, binwidth, smoothwidth)
