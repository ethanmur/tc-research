import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import eyewall_slope

import tc_metadata


# only two x axis choices: 'dist' and 'lat-lon' here! the xaxis touples above
# specify whether to use lat or lon
def plot_tdr_only( tc='sam', eyewall_cutoff=True, xaxis='dist'):

    tcdata = tc_metadata.choose_data_eyewall_slope( tc)
    if tcdata == 'selected TC name is not yet implemented':
        print( tcdata)
        return

    # print( 'number of tdr files: ' + str( len( tcdata['tdr_list'] )))
    for counter in range( len( tcdata[ 'dates'] )):
        plot_one_tdr_section( tcdata, counter, eyewall_cutoff, xaxis)

        # this code makes a blank figure in between the real figures: it basically
        # just makes the output look prettier
        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(.75) # figure height in inches


# a helper function to split up the plot_tdr_only() code!
# It also has some uses in the cloud_height_auto_same_plot.py file!
def plot_one_tdr_section( tcdata, counter, eyewall_cutoff, xaxis, good_data_case=False):
    warnings.filterwarnings("ignore")

    # load data
    os.chdir( tcdata['tdr_path'] )

    # special case for only using certain good datasets(good_data_case is set to true in the parent function)
    if good_data_case:
        inbound_data, outbound_data = tc_metadata.choose_tdr_data( tcdata['tc_name'], tcdata['tdr_list'], counter)
    # special case for ida (because most of the data looks so bad :(   )
    elif tcdata['tc_name'] == 'Ida':
        print( 'ida')
        inbound_data = tcdata[ 'tdr_list'] [ 20]
        outbound_data = tcdata['tdr_list'] [ 21]
    else:
        inbound_data = tcdata[ 'tdr_list'] [ counter*2]
        if len( tcdata['tdr_list']) > counter*2+1:
            outbound_data = tcdata['tdr_list'] [ counter*2 + 1]
        else:
            print( 'error!')

    if xaxis == 'dist':
        make_plots.plot_tdr( tcdata['tdr_path'], inbound_data, outbound_data, xaxis)
        x_in, x_out, H_in, H_out = eyewall_slope.eyewall_slope_first_val( tcdata['tdr_path'], inbound_data, outbound_data, eyewall_cutoff, xaxis)
        instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis)
    # lat case ( all lats are greater than 0 here)
    elif tcdata['xlims'][ counter][0] > 0:
        make_plots.plot_tdr( tcdata['tdr_path'], inbound_data, outbound_data, xaxis= 'lat')
        x_in, x_out, H_in, H_out = eyewall_slope.eyewall_slope_first_val( tcdata['tdr_path'], inbound_data, outbound_data, eyewall_cutoff, xaxis='lat')
        instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis='lat')
    # lon case
    else:
        x_in, x_out, H_in, H_out = eyewall_slope.eyewall_slope_first_val( tcdata['tdr_path'], inbound_data, outbound_data, eyewall_cutoff, xaxis='lon')
        instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis='lon')
        make_plots.plot_tdr( tcdata['tdr_path'], inbound_data, outbound_data, xaxis= 'lon')
    # the following comment is done in the code above:
    # add lines and red dots representing eyewall slope and beginning of eyewall, respectively
    # This follows the same procedure as the one used in TDR eyewall finder script.py!
    # calculate and load eyewall slope data


    # plot inbound eyewall data
    # first two calls are to plot eyewall slope until cutoff, third call plots point at start of eyewall
    plt.scatter(  x_in, H_in, c='b', s=20)
    plt.plot(  x_in, H_in, c='b', linewidth=1)

    # plt.scatter( - instartx, instarth, c='r', s=40) # older

    # plot a line marking the start of the eyewall at 50 km if no eyewall can be found
    if np.isnan( instartx):
        plt.axvline( x =  50, linewidth=2, c='g')
    else:
        plt.axvline( x= instartx, linewidth=2, c='g')

    # plot outbound eyewall data
    plt.scatter( x_out, H_out, c='b', s=20)
    plt.plot( x_out, H_out, c='b', linewidth=1)
    # plt.scatter( - outstartx, outstarth, c='r', s=40)
    if np.isnan( outstartx):
        plt.axvline( x = - 50, linewidth=2, c='g')
    else:
        plt.axvline( x=  outstartx, linewidth=2, c='g')

    if xaxis == 'dist':
        plt.xlim( [-100, 100])
    else:
        plt.xlim( [ tcdata['xlims'][ counter][0], tcdata['xlims'][counter][1] ] )
        # plt.xlim( [instartx, outstartx])

    os.chdir( "/Users/etmu9498/research/figures/tdr/eyewall-slopes/" + tcdata['tc_name'].casefold() + "/")
    if xaxis == 'dist':
        plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
    else:
        plt.savefig( tcdata['tc_name'].casefold() + "-lat-lon-" + str( counter+1) + ".png" )

    # print( "TDR Image " + str( counter + 1) + " complete" )

    warnings.filterwarnings("default")
