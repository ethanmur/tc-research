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
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here

    if tc.casefold() == 'fred':
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/fred/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        dates = [
            '08-11', '08-11',
            '08-12-am', '08-12-am', '08-12-am',
            '08-12-pm', '08-12-pm',
            '08-13', '08-13' ]
        eye_pass = ["1", "2", "1", "2", "3", "1", "2", "1", "2"]
        tc_name = 'Fred'

    elif tc.casefold() == 'grace':
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)


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

        dates = [
            '08-20', '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        eye_pass = ["1", "2", "3", "1", "2", "3"]
        tc_name = 'Henri'

    elif tc.casefold() == 'ida':
        # paths to data
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/ida/nc-files"
        # load a list of available data to these variables
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        dates = [
            '08-27', '08-27', '08-27', '08-27', '08-27', '08-27', '08-27',
            '08-28', '08-28', '08-28', '08-29', '08-29', '08-29']
        eye_pass = [
            "1", "2", "3","4","5","6","7",
            "1", "2", "3", "1", "2", "3"]
        tc_name = 'Ida'

    elif tc.casefold() == 'sam':
        # paths to data
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/sam/nc-files"
        # load a list of available data to these variables
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-27', '09-29', '09-29']
        eye_pass = [ "1", "2", "3", "1", "2", "3", "1", "2"]
        tc_name = 'Sam'

    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'crl_list': crl_list,
        'tdr_list': tdr_list, 'dates': dates, 'eye_pass': eye_pass, 'tc_name': tc_name }
    return tcdata




def plot_tdr_only( tc='sam', yaxis_zoom=False):

    warnings.filterwarnings("ignore")

    tcdata = choose_data( tc)
    if tcdata == 'selected TC name is not yet implemented':
        print( tcdata)
        return

    print( 'number of tdr files: ' + str( len( tcdata['tdr_list'] )))
    for counter in range( len( tcdata[ 'dates'] )):

        axis = 'distance'

        # load data
        os.chdir( tcdata['tdr_path'] )
        inbound_data = tcdata[ 'tdr_list'] [ counter*2]
        if len( tcdata['tdr_list']) > counter*2+1:
            outbound_data = tcdata['tdr_list'] [ counter*2 + 1]
        else:
            print( 'error!')

        plt.figure(figsize=(18,4), facecolor='w')
        plt.title( "TDR Cross Section, TC " + tcdata['tc_name'] + ", "
            + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])
        make_plots.plot_tdr( tcdata['tdr_path'], inbound_data, outbound_data, axis)

        # add lines and red dots representing eyewall slope and beginning of eyewall, respectively
        # This follows the same procedure as the one used in TDR eyewall finder script.py!
        # calculate and load eyewall slope data
        x_in, x_out, H_in, H_out = eyewall_slope.eyewall_slope_first_val( tcdata['tdr_path'], inbound_data, outbound_data)
        instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data)

        # plot inbound eyewall data
        # the negative signs are to flip data so they're in line correctly
        plt.scatter( - x_in, H_in, c='b', s=20)
        plt.plot( - x_in, H_in, c='b', linewidth=1)
        plt.scatter( - instartx, instarth, c='r', s=40)

        # plot outbound eyewall data
        plt.scatter( - x_out, H_out, c='b', s=20)
        plt.plot( - x_out, H_out, c='b', linewidth=1)
        plt.scatter( - outstartx, outstarth, c='r', s=40)

        plt.xlim( [-100, 100])

        os.chdir( "/Users/etmu9498/research/figures/tdr/eyewall-slopes/" + tcdata['tc_name'].casefold() + "/")
        plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
        print( "Image " + str( counter + 1) + " complete" )

        # this code makes a blank figure in between the real figures: it basically
        # just makes the output look prettier
        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(2) # figure height in inches

    warnings.filterwarnings("default")
