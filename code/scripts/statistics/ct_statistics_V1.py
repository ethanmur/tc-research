# this function is kinda outdated... look in cloud-height-pdfs.py for more
# current data analysis methods!

import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import eyewall_slope
import cloud_height
import eyewall_slope_auto
import tc_metadata


def tc_eyes( tc='all'):

    return


def one_tc_eye( crl_path, crl_name, tdr_path, inbound_name, outbound_name, xaxis, cutoff_power, i1, i2 ):

    in_x_start, out_x_start, in_H_start, out_H_start = eyewall_slope.eyewall_start( tdr_path, inbound_name, outbound_name, xaxis)

    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)
    if xaxis == 'lon':
        xaxis_data = crl_data.Lon.values[i1:i2]
    elif xaxis == 'lat':
        xaxis_data = crl_data.Lat.values[i1:i2]
    else:
        print( "update the if statement on line 37 in ct_statistics_v1.py!")
        return

    i1, i2 = find_new_inds( xaxis_data, in_x_start, out_x_start)

    H, xaxis_value = cloud_height.find_cloud_heights(crl_name, cutoff_power, i1, i2, xaxis = xaxis)

    # plot the crl figure for helpful testing
    plt.figure( figsize=(18, 5))
    make_plots.plot_power_ch1( crl_path, crl_name, i1, i2, xaxis)
    plt.scatter( xaxis_value, H, c= 'r', s=8, marker='s')
    plt.plot( xaxis_value, H, c=  'r', linewidth=1, label= 'Cloud Top Height')
    plt.xlabel( 'Time (UTC, Hours)')
    plt.title( "TC Sam, 09/26/22, Eye 1, CRL Backscattered Power and Cloud Top Heights")

    # print out some basic stats from the height dataset
    print( "Number of data points:  " + str( i2 - i1))
    print( "Height value range:     " + str( H.min().values) + " km to " + str( H.max().values) + " km")
    print( "Height value mean:      " + str( H.mean().values) + " km")
    print( "Height value median:    " + str( H.median().values) + " km")

    # view simple histogram of cloud top heights!

    # **** see radiative transfer class, lab 1 or something, for how to make a histogram with a curved line superimposed
    #      over the points! You need to use another library, not matplotlib ****

    nbins = 15 # number of bins

    plt.figure( figsize=( 10, 8))
    plt.hist( H, nbins, density=True, facecolor='c', histtype = 'bar') # 'bar', 'barstacked', 'step', 'stepfilled'
    plt.axvline( x=H.mean(), c='g', label="Mean height value", linewidth=3)

    plt.ylabel( 'Probability of a Given Height Occurence')
    plt.xlabel( 'Height from Surface (Km)')
    # plt.title( "TC Sam, 09/26/22, Eye 1, Cloud Top Height Histogram")
    plt.grid(True)
    plt.legend()


def find_new_inds( xaxis_data, in_x_start, out_x_start):
    # define two lambda functions to find the closest crl values to the tdr eyewall
    crl_lower_lim = min( xaxis_data, key= lambda list_value : abs(list_value - in_x_start.values))
    crl_upper_lim = min( xaxis_data, key= lambda list_value : abs(list_value - out_x_start.values))

    # find the indices for the closest matching values to xaxis
    crl_lower_lim_i = np.where( xaxis_data == crl_lower_lim)
    crl_upper_lim_i = np.where( xaxis_data == crl_upper_lim)

    # crl_lower_lim_i is annoyingly stored as a touple, simplify this and return
    # the results
    i1 = crl_lower_lim_i[0][0]
    i2 = crl_upper_lim_i[0][0]

    return i1, i2
