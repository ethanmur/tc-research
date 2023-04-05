# import...
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import warnings
import sys
import pandas as pd
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023 as helper_fns

def pull_all_datasets( tc='all'):

    # load fl names for processed datasets, using the typical dictionary structure
    fl_data_path = "/Users/etmu9498/research/data/in-situ-noaa-processed/"
    orig_data_path = "/Users/etmu9498/research/data/"

    # do this for all flight level datasets
    if tc == 'all':
        # make a list of years where flight level data is present
        os.chdir( orig_data_path)
        yearlist = [name for name in os.listdir('in-situ-noaa-processed')
            if os.path.isdir(os.path.join('in-situ-noaa-processed', name))]
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        for yeari, yearval in enumerate( yearlist):
            filelist.append( helper_fns.load_flight_level( fl_data_path + yearval, print_files=False) )
    # do this for just one year
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):
        yearlist = [ tc]
        filelist = [ helper_fns.load_flight_level( fl_data_path + tc, print_files=False)]
    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")

    # load the full ships dataset
    print( "Loading SHIPS dataset")
    os.chdir(  "/Users/etmu9498/research/data/ships/")
    file1 = open('lsdiaga_1982_2021_sat_ts_5day.txt', 'r')
    Lines = file1.readlines()
    # there are a bunch of lines in the total ships file!!
    print( "SHIPS dataset created")
    print( "Number of lines in dataset: " + str( len( Lines)))

    # pull ships data for all years and files given as input
    for yeari, year in enumerate( yearlist):
        fl_path = "/Users/etmu9498/research/data/in-situ-noaa-processed/" + year

        if year == ['2022']:
            print("SHIPS data for 2022 has not yet been created. Skipping these cases for now.")
            continue
        for namei, fl_name in enumerate( filelist[ yeari]):

            print( "flight level data accessed: " + fl_name)

            os.chdir( fl_path)
            data = xr.open_dataset( fl_name)

            # find the median time for this dataset
            median = data.time[0].values + ( data.time[-1].values - data.time[0].values) / 2
            #print( 'first time: ' + str( data.time[0].values))
            #print( 'last time: ' + str( data.time[-1].values))
            #print( 'midpoint: ' + str( median))

            # find the closest 6 hour SHIPS interval to this time
            hours = np.array( [0, 6, 12, 18, 24, 30]) # add 24 and 30 hours as wraparound cases
            closetime = hours[ np.argmin( np.abs( hours - median)) ]

            # roll back wraparound cases to the correct time (no 24 or 30 h cases in SHIPS database!)
            wraparound = False
            if closetime > 18:
                wraparound = True
                closetime -= 24

            # print out useful metadata for this case!
            strdate = fl_name[4:8] # original date: correct if wraparound time selected
            # get the formatting correct
            strtime = str( closetime)
            if strtime == '0' or strtime == '6':
                strtime = '0' + strtime
            strname = fl_name[ 11:15].upper()
            # cut off the string at the underscore! allow for 2 and 3 character tc names
            strtemp = ''
            for chari, charval in enumerate( strname):
                if charval == '_':
                    break
                else:
                    strtemp += charval
            strname = strtemp

            print( '\n')
            print( 'name: ' + strname)
            print( 'date: ' + strdate)
            print( 'time: ' + strtime)
            print( 'wraparound time? ' + str( wraparound))

            header_inds = []
            flind = 0
            # go through all the lines
            for ind in range( len( Lines)):
                # get the heading lines, and look for this TC's cases!
                if 'HEAD' and strname  in Lines[ ind]:
                    # only keep correct year cases
                    if Lines[ ind][ 6 : 8] == str( year[ 2:4]):
                        # print out the header and its index!
                        # print( 'index = ' + str( ind))
                        # print( Lines[ ind])
                        # save this index in a separate list! add all indices for this tc + year combo, useful for wraparound times
                        header_inds.append( ind)
                        # find the cases with the valid date and time!
                        if strdate + ' ' + strtime in Lines[ ind]:
                            # no matter the case, find the index!
                            flind = ind

            # no valid dataset found case... just skip to the next case!
            if flind == 0:
                print("No SHIPS dataset found for this case. Skipping to the next case.\n")
                continue

            # account for wraparound cases!
            # wraparound cases: a bit annoying. Skip ahead 24 hours to the next date!
            # this work
            if wraparound:
                # index in the limited header_inds list of the header of interest!
                small_list_ind = header_inds.index( flind)
                # shift 24 hours (6h * 4) ahead in indices
                flind = header_inds[ small_list_ind + 4]
            # print results!
            # print( header_inds)
            # print( flind)
            # print( Lines[ flind])
            # add the corrected date here!
            correct_date = Lines[ flind][ 6 : 12]


            # get 24h old and current vmax values
            vmaxold = 'nan'
            vmax = 'nan'
            # 24h before
            # make sure that there's SHIPS data 24h before the flight! if not, we can't get the old intensity
            # idea: maybe when there's no - 24h data, but there is - 12h data, use the -12h column in that dataset! SHIPS usually has
            # that data saved :)


            if header_inds.index( flind) - 4 >= 0:
                header_ind = header_inds[ header_inds.index( flind) - 4]
                # only look in this date's variables!
                i0 = header_ind
                i1 = header_inds[ header_inds.index( flind) - 3]
                # search for vmax!
                for i in range( i0,  i1 ):
                    if 'VMAX' in Lines[i]:
                        vmaxold = str( int( Lines[i][12 : 15])) # the last 3 vals for vmax at 0 hours
                        break
            else:
                print( "No valid vmax data 24h before flight...")

            # current max values
            i0 = flind
            i1 = header_inds[ header_inds.index( flind) + 1]
            for i in range( i0, i1):
                if 'VMAX' in Lines[i]:
                    vmax = str( int( Lines[i][12 : 15])) # the last 3 vals for vmax at 0 hours
                    break

            # print( vmaxold)
            # print( vmax)

            # pull previous and current shear dirs and mags
            dirold = 'nan'
            direc = 'nan'
            magold = 'nan'
            mag = 'nan'

            # 24h before
            if header_inds.index( flind) - 4 >= 0:
                i0 = header_inds[ header_inds.index( flind) - 4]
                i1 = header_inds[ header_inds.index( flind) - 3]
                # search for vmax!
                for i in range( i0,  i1 ):

                    #if i == i0 or i == i0 + 1:
                        #print( Lines[ i])

                    if 'SHRD' in Lines[i]:
                        #print( Lines[ i])
                        magold = str( np.round( float( Lines[i][11 : 15]) / 10, 2)) # the last 4 vals for shrd at 0 hours
                    if 'SHTD' in Lines[i]:
                        #print( Lines[ i])
                        dirold = str( int( Lines[i][11 : 15])) # the last 4 vals for shtd at 0 hours
            else:
                print( "No valid shear data 24h before flight...")

            # current time
            i0 = flind
            i1 = header_inds[ header_inds.index( flind) + 1]
            # search for vmax!
            for i in range( i0,  i1 ):
                #if i == i0 or i == i0 + 1:
                #    print( Lines[ i])

                if 'SHRD' in Lines[i]:
                    #print( Lines[ i])
                    # divide by 10 to get to kts, see SHIPS documentation
                    mag = str(  np.round( float( Lines[i][11 : 15]) / 10, 2)) # the last 4 vals for shrd at 0 hours
                if 'SHTD' in Lines[i]:
                    #print( Lines[ i])
                    direc = str( int( Lines[i][11 : 15])) # the last 4 vals for shtd at 0 hours


            # new new and previous vmax, shear, and shear dir variables to a text document!
            # get to the correct folder
            os.chdir("/Users/etmu9498/research/data/ships-output/" )
            # make sure the folder exists! if not, create it!
            output_folder = year
            if not os.path.isdir( output_folder):
                os.makedirs( output_folder)
                print( 'New folder created: ' + output_folder)
            os.chdir("/Users/etmu9498/research/data/ships-output/" + year)

            # create the file
            file1 = open( fl_name[ : - 13] + '.txt', 'w+')

            # add the metadata!
            addlines = [ "SHIPS Date:\n", correct_date + '\n', "SHIPS Time (Hours, UTC):\n", strtime + '\n', "24h old vmax (kt):\n",
                        vmaxold + '\n', "Current vmax (kt):\n", vmax + '\n', "24h old shear direction (degrees):\n", dirold + '\n',
                        'Current shear direction (degrees):\n', direc + '\n', '24h old shear magnitude (kt):\n', magold + '\n',
                        'Current shear magnitude (kt):\n', mag + '\n']
            file1.writelines( addlines)
            file1.close()
            print("Important SHIPS data saved under data/ships-output/")
