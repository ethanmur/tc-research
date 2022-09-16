import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib

def display_data_files( path, data_type, print_files):
    """
    This function simply prints out and returns a list of data files found at a
    specified path.
    :param path: the user specified path to data files. This argument must always
                 be present.
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


# this quick function updates the sizes of fonts in plots, and it decreases clutter in code!
def change_font_sizes(small=14, medium=18 ):
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize

# arg1 must be an xarray dataset, and arg2 can be a number or a dataset with matching
# axes lengths! arg2 is a little more flexible
# this returns two arguments: the index where the lowest value is found, and the lowest value itself!
def xr_closest_val( arg1, arg2):

    index = ( np.abs(arg1 - arg2)).argmin().values
    value = arg1[ index].values
    return index, value

# this code should work for numpy datasets and numbers in arg1!
def closest_val( arg1, arg2):

    index = ( np.abs( np.subtract( arg1, arg2))).argmin()
    value = arg1[ index]
    return index, value

def add_blank_colorbar():
    # add an empty colorbar to make everything fit in line... kinda a
    # messy solution but it's ok for now!
    viridis = cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    white = np.array([ 1, 1, 1, 1])
    newcolors[:, :] = white
    white_cmap = ListedColormap(newcolors)

    map = matplotlib.cm.ScalarMappable(cmap= white_cmap, norm=matplotlib.colors.Normalize( vmin= 0, vmax= 1))
    cbar = plt.colorbar( mappable= map)
    cbar.set_ticks([])
    cbar.outline.set_visible(False)

def white_background():
    plt.rcParams['axes.facecolor']='default'
    plt.rcParams['savefig.facecolor']='w'

def mpl_defaults():
    plt.rcdefaults()
