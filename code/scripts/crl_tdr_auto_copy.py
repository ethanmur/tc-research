## last edited: 6/23/22
## Status: Automatically plots TDR and CRL cross sections with the same lat or lon coords.
##         This is working much better! some comparisons look perfect, others seem slightly off still (see 09/26 eye 1)

import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots

# the auto plot function defaults to
# other valid storm selections include 'all', 'fred', 'grace', 'henri', and 'ida'
def plot( tc='sam'):

    warnings.filterwarnings("ignore")

    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here
    if tc.casefold() == 'sam':
        # paths to data
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/sam/nc-files"
        # load a list of available data to these variables
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        # list of touples representing the min and max indices used to clip crl data
        crl_range = [
            (0, 2000), (2350, 3800), (4400, 5600), (0, 2000), (2600, 3900),
            (4800, 6500), (0, 2000), (2100, 3800)]
        xlims = [
            (-52, -49), (-51.5, -49.75), (-51.75, -50), (-55, -51), (18, 15),
            (-54, -52), (-59.5, -56.5), (-59.5, -56.5) ]

        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-27', '09-29', '09-29']
        eye_pass = [ "1", "2", "3", "1", "2", "3", "1", "2"]
        tc_name = 'Sam'

    elif tc.casefold() == 'grace':
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        crl_range = [
            (3000, 4500), (5000, 6600), # 8/16
            (0, 1800), (2250, 4100), (4800, 6550), # 8/17
            (0, 2100), (2800, 4400), (5000, 6200), # 8/18
            (0, 1800), (2100, 3900), (4000, 5700) ] # 8/19
        xlims = [
            ( -73, -69), ( -72, -69 ),
            ( 20, 16), ( -79, -74), (20, 17),
            ( 22, 18), (-87, -82.5 ), (21, 18 ),
            (-92, -89), (-92, -90), (-93, -89) ]

        dates = [
            '08-16', '08-16',
            '08-17', '08-17', '08-17',
            '08-18', '08-18', '08-18',
            '08-19', '08-19', '08-19' ]
        eye_pass = ["1", "2", "1", "2", "3", "1", "2", "3", "1", "2", "3"]
        tc_name = 'Grace'

    elif tc.casefold() == 'henri':
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/henri/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        crl_range = [
            (400, 1800), (2600, 4300), (5000, 6100), # 8/20
            (0, 2400), (3100, 5100), (6700, 8100) ] # 8/21
        xlims = [
            ( 29.5, 33.5), (-76, -72 ), (-75, -72 ),
            (33.5, 39), (-73, -67.5), (-72, -67) ]

        dates = [
            '08-20', '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        eye_pass = ["1", "2", "3", "1", "2", "3"]
        tc_name = 'Henri'
    else:
        print( "not yet implemented")
        return

    print( 'number of tdr files: ' + str( len( tdr_list )))
    for counter in range( len( crl_range)):

        if xlims[ counter][0] > 0.0:
            axis = 'lat'
        else:
            axis = 'lon'

        # load data
        os.chdir( tdr_path)
        inbound_data = tdr_list[ counter*2]
        outbound_data = tdr_list[ counter*2 + 1]
        os.chdir( crl_path)

        if dates[ counter] == '08-16':
            crl_data = crl_list[ 4]
        elif dates[ counter] == '08-17':
            crl_data = crl_list[ 6]
        elif dates[ counter] == '08-18':
            crl_data = crl_list[ 7]
        elif dates[ counter] == '08-19':
            crl_data = crl_list[ 8]
        elif dates[ counter] == '08-20':
            crl_data = crl_list[ 9]
        elif dates[ counter] == '08-21':
            crl_data = crl_list[ 11]
        elif dates[ counter] == '09-26':
            crl_data = crl_list[ 16]
        elif dates[ counter] == '09-27':
            crl_data = crl_list[ 17]
        elif dates[ counter] == '09-29':
            crl_data = crl_list[ 18]
        else:
            print( 'update if statement!')
            return

        plt.figure(figsize=(18,12), facecolor='w')
        plt.subplot(311)
        plt.title( "TDR and CRL, Hurricane " + tc_name + ", " + dates[ counter] + ", Eye Pass " + eye_pass[ counter])
        make_plots.plot_tdr( tdr_path, inbound_data, outbound_data, axis)
        plt.xlim( xlims[ counter][0], xlims[counter][1])

        plt.subplot(312)
        make_plots.plot_power_ch1(crl_path, crl_data, crl_range[counter][0], crl_range[counter][1], axis)
        plt.xlim( xlims[ counter][0], xlims[counter][1])

        plt.subplot(313)
        make_plots.plot_T(crl_path, crl_data, crl_range[counter][0], crl_range[counter][1], axis)
        plt.xlim( xlims[ counter][0], xlims[counter][1])


        os.chdir( "/Users/etmu9498/research/figures/tdr/" + tc_name.casefold() + "/")
        plt.savefig( "tdr-crl-" + tc_name.casefold() + "-" + str( counter+1) + ".png" )
        print( "Image " + str( counter + 1) + " complete" )


        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(2) # figure height in inches

    warnings.filterwarnings("default")
