import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots
import eyewall_slope
import cloud_height
import eyewall_slope_auto

# tc='sam' means that the function defaults to running on sam's data if no tc
# name is provided!
def choose_data( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here

    if tc.casefold() == 'grace':
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
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/henri/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        crl_range = [
            (400, 1800), (2600, 4300), (5000, 6100), # 8/20
            (0, 2400), (3100, 5100), (6700, 8100) ] # 8/21
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
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/ida/nc-files"
        # load a list of available data to these variables
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        # list of touples representing the min and max indices used to clip crl data
        crl_range = [(1500, 4700)]
        xlims = [(-92, -87)]

        dates = ['08-29']
        eye_pass = [ "2"]
        tc_name = 'Ida'

    elif tc.casefold() == 'sam':
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

    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'crl_list': crl_list,
        'tdr_list': tdr_list, 'crl_range': crl_range, 'xlims': xlims,
        'dates': dates, 'eye_pass': eye_pass, 'tc_name': tc_name }
    return tcdata

# this function simply looks at the date provided in the tcdata dictionary and
# returns the appropriate name of the crl dataset
def choose_crl_date( date, crl_list):
    if date == '08-16':
        crl_data = crl_list[ 4]
    elif date == '08-17':
        crl_data = crl_list[ 6]
    elif date == '08-18':
        crl_data = crl_list[ 7]
    elif date == '08-19':
        crl_data = crl_list[ 8]
    elif date == '08-20':
        crl_data = crl_list[ 9]
    elif date == '08-21':
        crl_data = crl_list[ 11]
    elif date == '08-27':
        crl_data = crl_list[ 12]
    elif date == '08-28':
        crl_data = crl_list[ 13]
    elif date == '08-29':
        crl_data = crl_list[ 14]
    elif date == '09-26':
        crl_data = crl_list[ 16]
    elif date == '09-27':
        crl_data = crl_list[ 17]
    elif date == '09-29':
        crl_data = crl_list[ 18]
    else:
        print( 'update if statement!')
        return
    return crl_data


# the auto plot function defaults to
# other valid storm selections include 'all', 'fred', 'grace', 'henri', and 'ida'
def plot( tc='sam', tdr_crl=False):
    warnings.filterwarnings("ignore")

    tcdata = choose_data( tc)
    if tcdata == 'selected TC name is not yet implemented':
        print( tcdata)
        return

    print( 'number of crl files: ' + str( len( tcdata['xlims'] )))
    for counter in range( len( tcdata[ 'dates'] )):

        if tcdata[ 'xlims'] [ counter][0] > 0.0:
            axis = 'lat'
        else:
            axis = 'lon'

        # load data
        os.chdir( tcdata['tdr_path'] )
        # special case to plot the one crl tdr case that works for ida!
        if tc == 'ida-crl':
            inbound_data = tcdata[ 'tdr_list'] [ 20]
            outbound_data = tcdata['tdr_list'] [ 21]
        else:
            inbound_data = tcdata[ 'tdr_list'] [ counter*2]
            outbound_data = tcdata['tdr_list'] [ counter*2 + 1]
        os.chdir( tcdata['crl_path'] )

        crl_data = choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

        # code to find cloud heights and plot estimation of top heights
        cutoff_power = -30 # seems to be in the middle of decent power thresholds!
        i1 = tcdata['crl_range'][counter][0]
        i2 = tcdata['crl_range'][counter][1]
        H, time = cloud_height.find_cloud_heights( crl_data, cutoff_power, i1, i2, xaxis=axis)

        # call a function to find the start of the eyewall from the inbound and outbound datasets
        instartx, outstartx, instarth, outstarth = eyewall_slope.eyewall_start( tcdata['tdr_path'], inbound_data, outbound_data, xaxis= axis)

        # plot just crl data case
        if not tdr_crl:
            # make a figure and title
            plt.figure(figsize=(18, 5), facecolor='w')
            plt.title( "CRL Data, TC " + tcdata['tc_name'] + ", "
                + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

            # plot crl data
            make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)

            # plot cloud top height data
            plt.scatter( time, H, c= 'r', s=8, marker='o') # s
            plt.plot( time, H, c= 'r', linewidth=1.5, label= 'Cloud Top Height')
            leg = plt.legend()
            for line in leg.get_lines():
                line.set_linewidth(4.0)

            # set crl data limits based on eyewall_start() function above
            if not np.isnan( instartx) and not np.isnan( outstartx):
                plt.xlim( instartx, outstartx)
            # the original way to set limits: the ones used in choose_data() function!
            # plt.xlim( tcdata['xlims'] [counter][0], tcdata['xlims'] [counter][1])

            os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-auto/" + tcdata['tc_name'].casefold() + "/")
            plt.savefig( "crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
            print( "Image " + str( counter + 1) + " complete" )

        # plotting tdr and crl data case. Note: saves figs to a different location!
        else:
            # make a figure and title
            plt.figure(figsize=(18, 10), facecolor='w')
            plt.subplot(211)
            plt.title( "TDR and CRL, TC " + tcdata['tc_name'] + ", "
                + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])

            # plot tdr data
            eyewall_slope_auto.plot_one_tdr_section( tcdata, counter, eyewall_cutoff=False, xaxis = axis)

            step = .2

            if not np.isnan( instartx) and not np.isnan( outstartx):
                print( instartx.values)
                print( str(outstartx.values) + '\n' )
                print( 'tdr no nan')
                if instartx < outstartx:
                    plt.xlim( instartx - step, outstartx + step)
                else:
                    plt.xlim( instartx + step, outstartx - step)

            plt.subplot( 212)
            # plot crl data
            make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)


            # plot cloud top height data
            plt.scatter( time, H, c= 'r', s=8, marker='o') # s
            plt.plot( time, H, c= 'r', linewidth=1.5, label= 'Cloud Top Height')
            leg = plt.legend()
            for line in leg.get_lines():
                line.set_linewidth(4.0)

            # plot eyewall cutoffs
            plt.axvline( x= instartx, linewidth=2, c='g')
            plt.axvline( x= outstartx, linewidth=2, c='g')

            # set crl data limits based on eyewall_start() function above
            if not np.isnan( instartx) and not np.isnan( outstartx):
                print( 'crl no nan')
                if instartx < outstartx:
                    plt.xlim( instartx - step, outstartx + step)
                else:
                    plt.xlim( instartx + step, outstartx - step)

            # the original way to set limits: the ones used in choose_data() function!
            # plt.xlim( tcdata['xlims'] [counter][0], tcdata['xlims'] [counter][1])

            os.chdir( "/Users/etmu9498/research/figures/CRL-eye-profiles-auto/" + tcdata['tc_name'].casefold() + "/")
            plt.savefig( "tdr-crl-cloud-height-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
            print( "CRL Image " + str( counter + 1) + " complete" )

        # this code makes a blank figure in between the real figures: it basically
        # just makes the output look prettier
        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(2) # figure height in inches

    warnings.filterwarnings("default")
