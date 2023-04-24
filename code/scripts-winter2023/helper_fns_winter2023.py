import os
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors

# get and return a list of tc years and file names for a specific subset (or all) crl tc data!
# based off code from plotting in "code/2022 eyewall determinations.ipynb"
def get_crl_datasets( tc='all', crl_data_root = "/Users/etmu9498/research/data/crl-all-data-processed/"):

    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( load_flight_level( crl_data_root + '2022', print_files=False) )

    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ load_flight_level( crl_data_root + '2022', print_files=False)]

    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")

    return yearlist, filelist

# the same process as above, but for flight level data
def get_fl_datasets( tc='all', crl_data_root = "/Users/etmu9498/research/data/in-situ-noaa-processed/"):
    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']
        # make a list of lists of all the datasets to be processed ( each year has a sublist)
        filelist = []
        filelist.append( load_flight_level( crl_data_root + '2021', print_files=False) )
        filelist.append( load_flight_level( crl_data_root + '2022', print_files=False) )

    # do this for just one year
    elif tc == '2021':
        yearlist = ['2021']
        filelist = [ load_flight_level( crl_data_root + '2021', print_files=False)]
    elif tc == '2022':
        yearlist = ['2022']
        filelist = [ load_flight_level( crl_data_root + '2022', print_files=False)]

    # do this for a specific dictionary of files:
    elif type( tc) == type( {}):
        # make the folder_list: just the dates saved in the dictionary!
        yearlist = list( tc.keys())
        filelist = []
        for keyi, keyval in enumerate( yearlist):
            filelist.append( tc[ keyval])
    else:
        print( "Error: please enter a valid selection for tc")
    return yearlist, filelist


# similar to the process above, but for tdr data. also return a list of date folders ( a bit more complicated)
def get_tdr_datasets( tc='all', tdr_data_root = "/Users/etmu9498/research/data/tdr-original/"):
    # do this for all crl datasets (2021 and 2022)
    if tc == 'all':
        # make a list of years where crl data is present
        yearlist = ['2021', '2022']
        
        namelist = []
        sublist = []
        path = tdr_data_root + "2021/"
        for (dirpath, dirnames, file) in os.walk( path):
            sublist.extend(dirnames)
            break
        namelist.append( sublist)

        # repeat for 2022
        sublist = []
        path = tdr_data_root + "2022/"
        for (dirpath, dirnames, file) in os.walk( path):
            sublist.extend(dirnames)
            break
        namelist.append( sublist)


        # finally, create the file list of cases to look at
        filelist = [] 
        # empty lists for 2021, 2022
        # add more empty lists for each tc name and year
        for yeari, yearval in enumerate(yearlist):
            filelist.append([])
            for tci, tcval in enumerate(namelist[yeari]):
                filelist[yeari].append( [])
        print(filelist)

        for yeari, yearval in enumerate(yearlist):

            for filei, filename in enumerate(namelist[ yeari]):
                sublist = []
                path = tdr_data_root + yearval + "/" + filename + '/'
                for (dirpath, dirnames, file) in os.walk( path):
                    sublist.extend(file)
                    break
                filelist[yeari][filei].append(sublist)

    # do this for just one year
    # elif tc == '2021':
    #     yearlist = ['2021']
    #     filelist = [ load_flight_level( tdr_data_root + '2021', print_files=False)]
    # elif tc == '2022':
    #     yearlist = ['2022']
    #     filelist = [ load_flight_level( tdr_data_root + '2022', print_files=False)]
    else:
        print( "Error: please enter a valid selection for tc")

    return yearlist, namelist, filelist



# find tc intensity category
def find_category( intensity_list):
    cat = []
    for i in range( len( intensity_list)):
        # tropical depression case (all values are in knots)
        if intensity_list[ i] < 34:
        # tropical storm
            cat += ['td']
        elif intensity_list[ i] >= 34 and intensity_list[ i] < 64:
            cat += ['ts']
        # weak hurricane (categories 1-2)
        elif intensity_list[ i] >= 64 and intensity_list[ i] < 96:
            cat += ['wh']
        # strong hurricane (categories 3-5)
        elif intensity_list[ i] >= 96:
            cat += ['sh']
        else:
            print( 'Fix if statement in intensity_cat() function!')
    return cat

def change_font_sizes(small=14, medium=18 ):
    """this quick function updates the sizes of fonts in plots, decreasing code clutter.
    it only applies the font changes if called before figure / axis creation.

    param small: change font sizes to this smaller value.
    param medium: change font sizes to this larger value.
    return: none
    """
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize


def load_crl( path, print_files=True):
    """load crl data at the given data path using a helper function."""
    return display_data_files( path, 'crl', print_files)

def load_flight_level( path, print_files=True):
    """load flight level data at the given data path using a helper function."""
    return display_data_files( path, 'in-situ', print_files)

def display_data_files( path, data_type, print_files):
    """This function simply prints out and returns a list of data files found at a
    specified path.
    :param path: the user specified path to data files. This argument must always be present.
    :param data_type: must include a data identifier, like 'crl', 'tdr', 'in-situ',
    or 'dropsonde'. Can also include 'hide-list' to surpress printing
    a list of the data
    :param print_status: Implicitly print the file names if given 'show-list',
    unless "hide-list" is passed instead.
    :return: A list of file names found within the specified folder.
    """
    file_names = []
    for (dirpath, dirnames, file) in os.walk( path):
        file_names.extend(file)
        break

    if data_type == 'crl' and print_files:
        print( "crl data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'tdr' and print_files:
        print( "tdr data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'in-situ' and print_files:
        print( "in situ data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'dropsonde' and print_files:
        print( "dropsonde data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'goes' and print_files:
        print( "GOES satellite data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])

    return file_names



def add_blank_colorbar(color=False):
    # add an empty colorbar to make everything fit in line... kinda a
    # messy solution but it's ok for now!
    viridis = cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))

    if color == 'gray':
        white = np.array([ .92, .92, .92, .92])
    else:
        white = np.array([ 1, 1, 1, 1])

    newcolors[:, :] = white
    white_cmap = ListedColormap(newcolors)

    map = matplotlib.cm.ScalarMappable(cmap= white_cmap, norm=matplotlib.colors.Normalize( vmin= 0, vmax= 1) )
    cbar = plt.colorbar( mappable= map)
    cbar.set_ticks([])
    cbar.outline.set_visible(False)