import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import eyewall_slope

def choose_data( tc='sam'):

    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    crl_list = make_plots.load_crl(crl_path, print_files=False)
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here

    if tc.casefold() == 'fred':
        tdr_path = "/Users/etmu9498/research/data/tdr/fred/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        xlims = [
            ( -82.5, -68), ( -82.5, -68 ),
            ( -77, -72), ( -77, -72), (24, 18),
            ( -77, -73.5), (-77, -73.5 ),
            (-79, -75), (-79, -75) ]
        dates = [
            '08-11', '08-11',
            '08-12-am', '08-12-am', '08-12-am',
            '08-12-pm', '08-12-pm',
            '08-13', '08-13' ]
        eye_pass = ["1", "2", "1", "2", "3", "1", "2", "1", "2"]
        tc_name = 'Fred'

    elif tc.casefold() == 'grace':
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        xlims = [
            ( -73, -69), ( -72, -69 ),
            ( 20, 16), ( -79, -74), (20, 17),
            ( 22, 18), (-86.5, -82.5 ), (21, 18 ),
            (-92, -89), (-92, -90), (-93, -89) ]
        dates = [
            '08-16', '08-16',
            '08-17', '08-17', '08-17',
            '08-18', '08-18', '08-18',
            '08-19', '08-19', '08-19' ]
        eye_pass = ["1", "2", "1", "2", "3", "1", "2", "3", "1", "2", "3"]
        tc_name = 'Grace'

    elif tc.casefold() == 'henri':
        tdr_path = "/Users/etmu9498/research/data/tdr/henri/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        xlims = [
            ( 29.5, 33.5), (-76, -72 ), (-75, -72 ),
            (33.5, 39), (-73, -67.5), (-72, -68.5) ]
        dates = [
            '08-20', '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        eye_pass = ["1", "2", "3", "1", "2", "3"]
        tc_name = 'Henri'

    elif tc.casefold() == 'ida':
        # paths to data
        tdr_path = "/Users/etmu9498/research/data/tdr/ida/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        xlims = [
            (-85, -81), (-85, -81), (-85, -83), (24, 23), (-85, -83), (-85, -83),
            (-86.5, -82), (-88.5, -85.5), (-88, -86),
            (-90, -85), (-92, -87.5), (-92, -88.5), (-92, -88.5)]
        dates = [
            '08-27', '08-27', '08-27', '08-27', '08-27', '08-27', '08-27',
            '08-28', '08-28', '08-28', '08-29', '08-29', '08-29']
        eye_pass = [
            "1", "2", "3","4","5","6","7",
            "1", "2", "3", "1", "2", "3"]
        tc_name = 'Ida'

    elif tc.casefold() == 'sam':
        tdr_path = "/Users/etmu9498/research/data/tdr/sam/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        xlims = [
            (-52, -49), (-51.5, -49.75), (-51.75, -50), (-55, -51), (18, 15),
            (-54, -52), (-59.5, -56.5), (-59.5, -56.5) ]
        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-27', '09-29', '09-29']
        eye_pass = [ "1", "2", "3", "1", "2", "3", "1", "2"]
        tc_name = 'Sam'
    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'crl_list': crl_list,
        'tdr_list': tdr_list, 'xlims': xlims, 'dates': dates,
        'eye_pass': eye_pass, 'tc_name': tc_name }
    return tcdata



# only two x axis choices: 'dist' and 'lat-lon' here! the xaxis touples above
# specify whether to use lat or lon
def plot_tdr_only( tc='sam', eyewall_cutoff=True, xaxis='dist'):

    tcdata = choose_data( tc)
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
def plot_one_tdr_section( tcdata, counter, eyewall_cutoff, xaxis):
    warnings.filterwarnings("ignore")

    # load data
    os.chdir( tcdata['tdr_path'] )
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
    print( "TDR Image " + str( counter + 1) + " complete" )

    warnings.filterwarnings("default")
