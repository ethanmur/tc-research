import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

import warnings
import pandas as pd
from scipy.signal import find_peaks

os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots_new_heights
import make_plots
import tc_metadata
import helper_fns

os.chdir("/Users/etmu9498/research/code/scripts/in-situ-scripts")
import load_in_situ_data

import sys
sys.path.append("/Users/etmu9498/research/code/scripts/plotting")
import rmw_inner_outer_stats as extra_stats # confidence interval function resides here!
import fl_vpeaks_algorithms as peak_algs


# this function plots all eye passes from RMW = -1 to RMW = 1
# fl data hasn't been binned at this point- more just a test to make sure data looks good
def plot_all_eyes(tc='all', max_v_requirement=40, window=10, method='pmin_tlim', timelim=30, filter_inner_core=False):
    # load data
    # this code comes from auto_flight_level_plots_new_noaa_data.py
    fl_path_root = "/Users/etmu9498/research/data/in-situ-noaa-full/"
    data_path = "/Users/etmu9498/research/data/"

    # make a list of year folders saved in 'in-situ-noaa-full'
    os.chdir( data_path)
    folder_list = [name for name in os.listdir('in-situ-noaa-full')
        if os.path.isdir(os.path.join('in-situ-noaa-full', name))]
    fl_list_total = []
    # go through every year- count the number of datasets!
    for folderi in range( len( folder_list)):
        # go to new directory and get the count + names of files there!
        fl_listi = make_plots.load_flight_level( fl_path_root + folder_list[ folderi], print_files=False)
        fl_list_total += fl_listi
    fl_list_count = len( fl_list_total)


    # use the whole list, or make plots for just one year!
    if tc == 'all':
        print( 'Total Number of plots to be created: ' + str(fl_list_count)+ '\n')
        folder_input = folder_list # use all possible year folders

    # one year case- check the first two characters: a proper year input?
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):
        # check if this year has a folder
        folder_present = False
        for folder in folder_list:
            # the folder exists
            if folder == tc:
                folder_input = [ tc]
                folder_present = True
                fl_listi = make_plots.load_flight_level( fl_path_root + tc, print_files=False)
                print( 'Total Number of plots to be created: ' + str( len( fl_listi))+ '\n')
        # no folder present case:
        if folder_present == False:
            print( "No folder present. Please download flight level data and add a new folder for it!")
            return

    # input names of specific runs as a dict! each year has its own entry
    elif type( tc) == type( {}):

        # make the folder_list: just the dates saved in the dictionary!
        folder_input = list( tc.keys())

    # made it out of the loop?
    else:
        print( 'implement the else case! cannot handle individual plots yet.')




    # initialize the dataframe!!
    fl_df_total = pd.DataFrame( )
    # do this for every year applicable: all years, one year, or a special case!
    for yeari in folder_input:

        # make a temp dataframe that will be concatonated to the final dataframe!
        fl_df_current = pd.DataFrame( )

        print( "\nCurrent year: " + yeari)


        fl_path = fl_path_root + yeari

        # normal cases: just load all the datasets from the given year!
        if tc=='all' or ( len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20')):
            fl_list = make_plots.load_flight_level( fl_path, print_files=False)
        # dictionary case: load only the datasets provided!
        elif type( tc) == type( {}):
            fl_list = tc[ yeari]

        # new code
        # add the automatic variable names to the dataset
        auto_vars = [ 'name', 'year', 'pass', 'rmw', 'Time'] # always find these variables
        # input variables
        input_vars = [ 'WS.d', 'UWZ.d', 'SfmrRainRate.1', 'THETAE.d', 'MR.d', 'TA.d'] # change me!
        var_names = [ 'wind_speed', 'w', 'Rain Rate', 'Theta E', 'Mixing Ratio', 'temp']  # change me!





        # load datasets once here to save time / repetition :)
        os.chdir( fl_path)
        data_list = []

        fl_list_iterate = []
        for i in range( len( fl_list)):
            # load data
            dataset = xr.open_dataset( fl_list[ i], decode_times=False)
            dataset = process_data_step1( dataset)

            # smooth wind speeds and surface pressures- more realistic peaks
            window = window # seconds
            spd_avg = pd.Series( dataset['WS.d']).rolling(window=window, min_periods=1, center=True).mean()
            # test 1: wind speed threshold with no eyewall limits
            # find max TC strength from averaged wind fields- all eye passes from this tc
            vmax = np.nanmax( spd_avg)
            # check based on max TC strength
            if vmax > max_v_requirement:
                data_list.append( dataset)
                fl_list_iterate.append( fl_list[ i])
                print( "Dataset " + str( i) + " loaded")

            else:
                print( "Case " + str( i ) + " not included: peak winds < " + str( np.round( max_v_requirement, 1)) + " m/s")

        fl_list = fl_list_iterate
        #print( fl_list)

        #############################
        # new code
        # add the fixed / automatic variables to the dataframe first
        #############################
        peak_min = 40


        for input_ind in range( len( auto_vars)):
            var_total = []

            input_field = auto_vars[ input_ind]

            # print( 'Variable added: ' + input_field)

            # loop through flight days
            for i in range( len( fl_list)):
                dataset = data_list[ i]

                # smooth wind speeds and surface pressures- more realistic peaks
                window = window # seconds
                spd_avg = pd.Series( dataset['WS.d']).rolling(window=window, min_periods=1, center=True).mean()

                # find max vel peaks using helper functions
                if method == 'pmin':
                    vpeaks, pmins = peak_algs.find_peaks_pressure_mins( dataset['PSURF.d'], spd_avg, window)
                elif method == 'pmin_tlim':
                    vpeaks, pmins, time_lims = peak_algs.find_peaks_pmin_time_limit( dataset['PSURF.d'], spd_avg,
                            dataset['Time'].values, window, timelim=timelim, filter_inner_core = filter_inner_core)

                pass_counter = 0
                # loop through eyewall passes
                # do step = 2 to skip the second index for each eyewall!
                for vind in range( int( len(vpeaks) / 2 ) ):

                    ###########################
                    ## new code
                    ## do different things depending on the input var!
                    ###########################
                    # quick cases:
                    if input_field == 'name':
                        var_total.append( fl_list[ i])
                    elif input_field == 'year':
                        var_total.append( fl_list[ i][0:4])
                    elif input_field == 'pass':
                        var_total.append( pass_counter)

                    # longer cases
                    elif input_field == 'Time':
                        left_eyewall_ind = vpeaks[ 2* vind]
                        right_eyewall_ind = vpeaks[2 * vind + 1]

                        # find the corresponding times
                        ws_ind = np.nanargmin( spd_avg[ left_eyewall_ind:right_eyewall_ind])

                        # trim down flight level data for a given field to just the rmw region and append it to the temp variable
                        var_total.append( dataset[ auto_vars[ input_ind]] [ left_eyewall_ind : right_eyewall_ind])

                    elif input_field == 'rmw':
                        left_eyewall_ind = vpeaks[ 2* vind]
                        right_eyewall_ind = vpeaks[2 * vind + 1]

                        # find the corresponding times
                        ws_ind = np.nanargmin( spd_avg[ left_eyewall_ind:right_eyewall_ind])

                        # make an rmw axis stretching from one eyewall to the other!
                        # they will have an unequal number of data points to conserve the time axis :)
                        len1 = ws_ind
                        len2 = right_eyewall_ind - ( left_eyewall_ind + ws_ind)
                        var_total.append( np.concatenate( [ np.linspace( -1, 0, num=len1, endpoint=False), np.linspace( 0, 1, num=len2) ] ))
                    pass_counter += 1
            # add this current variable to temp df (the indent is intentional!)
            fl_df_current[ auto_vars[ input_ind]] = var_total
        print( 'Auto variables added')







        ####################
        # new code
        # repeat this step for every chosen variable in the input list!
        # add a new column to the dataframe for each input variable
        ####################
        for input_ind in range( len( input_vars)):

            # total variable that saves individual runs for each tc
            # if there are 3 eye passes for 2 storms, var_total = [ [ temps], [temps2], [temps3] ]
            var_total = []

            input_field = input_vars[ input_ind]

            # print( 'Variable added: ' + input_field)

            # loop through flight days
            for i in range( len( fl_list)):
                dataset = data_list[ i]

                # smooth wind speeds and surface pressures- more realistic peaks
                window = window # seconds
                spd_avg = pd.Series( dataset['WS.d']).rolling(window=window, min_periods=1, center=True).mean()

                # test 1: pressure threshold with no eyewall limits
                # find max TC strength from averaged wind fields- all eye passes from this tc
                vmax = np.nanmax( spd_avg)

                # find max vel peaks using helper functions
                if method == 'pmin':
                    vpeaks, pmins = peak_algs.find_peaks_pressure_mins( dataset['PSURF.d'], spd_avg, window)
                elif method == 'pmin_tlim':
                    vpeaks, pmins, time_lims = peak_algs.find_peaks_pmin_time_limit( dataset['PSURF.d'], spd_avg, dataset['Time'].values, window, timelim=timelim)

                pass_counter = 0
                # loop through eyewall passes
                # do step = 2 to skip the second index for each eyewall!
                for vind in range( int( len(vpeaks) / 2 ) ):


                    left_eyewall_ind = vpeaks[ 2* vind]
                    right_eyewall_ind = vpeaks[2 * vind + 1]

                    # find the corresponding times
                    ws_ind = np.nanargmin( spd_avg[ left_eyewall_ind:right_eyewall_ind])

                    # make an rmw axis stretching from one eyewall to the other!
                    # they will have an unequal number of data points to conserve the time axis :)
                    len1 = ws_ind
                    len2 = right_eyewall_ind - ( left_eyewall_ind + ws_ind)

                    # trim down flight level data for a given field to just the rmw region
                    # and append it to the temp variable
                    var_total.append( dataset[ input_vars[ input_ind]] [ left_eyewall_ind : right_eyewall_ind])
                    pass_counter += 1

            # add new column to the dataframe!
            fl_df_current[ var_names[ input_ind]] = var_total

        print("Chosen variables added")
        # once finished making the current dataframe, concatonate it onto the total!
        fl_df_total = pd.concat( [fl_df_total, fl_df_current], ignore_index=True)


    # print( fl_df_total)
    fl_data_totals = fl_df_total
    ######################
    ## Make a nice rmw plot for each eye pass!
    ######################

    # make figure before loop and add to it!
    plt.figure( figsize=(14, 28))
    lw = 2.5
    helper_fns.change_font_sizes(16, 16)

    plt.subplot(611)
    plt.title( "All fields of interest from TC Ida and Sam")
    plt.ylabel("Total Wind Speed (m/s)")
    plt.xlim( [-1.2, 1.2])
    plt.grid()
    plt.subplot(612)
    plt.ylabel("Vertical Velocity (m/s)")
    plt.xlim( [-1.2, 1.2])
    plt.ylim( [-15, 35])
    plt.grid()

    plt.subplot( 613)
    plt.ylabel("Temp ( C)")
    plt.xlim( [-1.2, 1.2])

    plt.subplot( 614)
    plt.ylabel("Theta E ( Kelvin)")
    plt.xlim( [-1.2, 1.2])
    plt.subplot(615)
    plt.ylabel("WVMR (g/kg))")
    plt.xlim( [-1.2, 1.2])
    plt.subplot(616)
    plt.ylabel("SFMR Rain Rate (mm/hr)")
    plt.xlabel("Radius of Maximum Winds (RMW)")
    plt.xlim( [-1.2, 1.2])


    # create a colormap for nice plotting!
    start = 0.0
    stop = 1.0
    number_of_lines= fl_data_totals.shape[0]

    cm_subsection = np.linspace(start, stop, number_of_lines)
    colors = [ cm.jet(x) for x in cm_subsection ]


    # do this for every flight date in this case (defined in code block above)
    # print( "Number of TC flight days = " + str( len( fl_list_total)))
    print( 'Total number of passes = ' + str( number_of_lines))

    # cycle through each row in the dataframe
    for i in range( number_of_lines):

        # make plots of total wind speed and vertical vels along new rmw axes!
        plt.subplot(611)
        plt.plot( fl_data_totals['rmw'][i], fl_data_totals['wind_speed'][i], c=colors[ i], linewidth=lw, alpha=.7)

        plt.subplot(612)
        plt.plot( fl_data_totals['rmw'][i], fl_data_totals['w'][i], c=colors[ i], linewidth=lw, alpha=.7)

        plt.subplot(613)
        label = "Flight " + fl_data_totals['name'][i] + ", pass " + str( fl_data_totals['pass'][i])
        plt.plot( fl_data_totals['rmw'][i], fl_data_totals['temp'][i], c=colors[ i], label=label, linewidth=lw, alpha=.7)
        if i == number_of_lines - 1:
            plt.legend( bbox_to_anchor=(1.03, 1.05), fancybox=False, shadow=False, fontsize=14, facecolor='w', framealpha=1)
            plt.grid()


        plt.subplot(614)
        plt.plot( fl_data_totals['rmw'][i], fl_data_totals['Theta E'][i], c=colors[ i], label=label, linewidth=lw, alpha=.7)

        plt.subplot(615)
        # trim out non physical peaks in wv from plots! maybe do this when data is saved?
        # this actually doesn't help, so I commented it out :(

        rmw_trim = fl_data_totals['rmw'][i] # [ np.where( fl_data_totals['Mixing Ratio'][i] < 30 ) ]
        wv_trim = fl_data_totals['Mixing Ratio'][i] # [ np.where( fl_data_totals['Mixing Ratio'][i] < 30 )]
        plt.plot( rmw_trim, wv_trim, c=colors[ i], linewidth=lw, alpha=.7)

        plt.subplot(616)
        plt.plot( fl_data_totals['rmw'][i], fl_data_totals['Rain Rate'][i], c=colors[ i], linewidth=lw, alpha=.7)






def process_data_step1( dataset):
    # this snippet turns all relevant quantities from strings into floats
    # add more fields here if needed!
    keylist = [ 'Time', 'WS.d', 'WD.d', 'UWZ.d', 'SfmrRainRate.1', 'LATref', 'LONref', 'TAS.d', 'PSURF.d',
               'HT.d', 'PITCHref', 'ROLLref', 'MR.d', 'HUM_REL.d', 'THETA.d', 'THETAE.d', 'THETAV.d', 'TA.d']
    # create a new version of the xarray dataset holding only the provided keys
    datasetnew = dataset.copy()
    datasetnew = datasetnew[ keylist]

    # create a new time axis!
    datasetnew[ 'time_index'] = dataset[ 'Time']

    # this string holds the start and end times. cut down the interval dataset to get the hours, mins, secs!
    interval_str = dataset.attrs['TimeInterval']
    h = float( interval_str[0:2])
    m = float( interval_str[3:5])
    s = float( interval_str[6:8])
    start_time = h + m / 60 + s / 3600

    # create the time array manually
    # if this is too slow, try using pandas? see stack overflow link below...
    # https://stackoverflow.com/questions/55648630/how-to-decode-the-time-variable-while-using-xarray-to-load-a-netcdf-file
    time_array = np.empty( ( len( dataset['Time'])))
    for timei in range( len( dataset['Time'])):
        # add to time array
        time_array[ timei] = start_time + timei / 3600
    datasetnew[ 'Time'] = time_array

    return datasetnew
