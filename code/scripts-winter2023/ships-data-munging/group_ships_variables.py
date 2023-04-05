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
os.chdir(  "/Users/etmu9498/research/code/scripts-winter2023/ships-data-munging")
import ships_plots

# sort through the text files and pull out important information! Like intensity,
# intensification categories.
def group_all( tc='all', make_histograms=False, ri_cond=30):
    # load ships names for processed datasets, using the typical dictionary structure
    ships_data_path = "/Users/etmu9498/research/data/ships-output/"
    orig_data_path = "/Users/etmu9498/research/data/"

    # do this for all flight level datasets
    if tc == 'all':
        # make a list of years where flight level data is present
        os.chdir( orig_data_path)
        yearlist = [name for name in os.listdir('ships-output')
            if os.path.isdir(os.path.join('ships-output', name))]
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        for yeari, yearval in enumerate( yearlist):
            filelist.append( helper_fns.load_flight_level( ships_data_path + yearval, print_files=False) )
    # do this for just one year
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):
        yearlist = [ tc]
        filelist = [ helper_fns.load_flight_level( ships_data_path + tc, print_files=False)]
    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")

    # store all shear, etc vals in one dataframe for easy acess!
    cols = [ 'Name', 'SHIPS Date', 'SHIPS Time', 'Vmax (kt)', 'Intensification (kt)', 'Shear (kt)', 'Shear Change (kt)']
    ships_df = pd.DataFrame( columns = cols)

    # pull ships data for all years and files given as input
    for yeari, year in enumerate( yearlist):
        fl_path = "/Users/etmu9498/research/data/ships-output/" + year
        os.chdir( fl_path)
        for namei, nameval in enumerate( filelist[ yeari]):
            # save local wind speeds, etc here
            vmax = np.nan
            shearmag = np.nan
            vchange = np.nan
            shearchange = np.nan

            # load the current text file
            file1 = open( nameval, 'r')
            Lines = file1.readlines()
            # go through all the text file lines
            for linei, lineval in enumerate( Lines):

                # find the date and time
                if lineval == "SHIPS Date:\n":
                    datei = Lines[ linei + 1][:-1]
                elif lineval == "SHIPS Time (Hours, UTC):\n":
                    timei = int( Lines[ linei + 1][:-1] )
                # find vmax info
                elif lineval == "Current vmax (kt):\n":
                    if Lines[ linei + 1][:-1] != 'nan':
                        vmax =  int( Lines[ linei + 1][:-1])# :-1 is to remove the \n from the string!
                    if Lines[ linei - 1][:-1] != 'nan' and Lines[ linei + 1][:-1] != 'nan':
                        vchange = int( Lines[ linei + 1][:-1]) - int( Lines[ linei - 1][:-1])
                # find shear info
                elif lineval == "Current shear magnitude (kt):\n":
                    if Lines[ linei + 1][:-1] != 'nan':
                        shearmag =  np.round( float( Lines[ linei + 1][:-1]), 3)# :-1 is to remove the \n from the string!
                    if Lines[ linei - 1][:-1] != 'nan' and Lines[ linei + 1][:-1] != 'nan':
                        shearchange = np.round( float( Lines[ linei + 1][:-1]) - float( Lines[ linei - 1][:-1]), 3)

            # add new variables to the dataframe!
            #ships_df.append( { 'Name' : nameval, 'Vmax (kt)' : vmax, 'Intensification (kt)' : vchange,
            #        'Shear (kt)' : shearmag, 'Shear Change (kt)' : shearchange}, ingnore_index=True)
            ships_df.loc[ len( ships_df.index)] = [ nameval, datei, timei, vmax, vchange, shearmag, shearchange]

    # optional code: print out pretty histograms of the intensities, intensity changes, etc!
    if make_histograms:
        ships_plots.plot_hists( ships_df, ri_cond)

    # return the final dataframe!
    return ships_df
