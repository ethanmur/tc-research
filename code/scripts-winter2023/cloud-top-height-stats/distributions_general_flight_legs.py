# import...
import numpy as np
import os
import sys
import pandas as pd
from geopy import distance
import metpy.calc as mpcalc
import xarray as xr
import math
import matplotlib.pyplot as plt
os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023


# Make a nice shear corrected figure of flight paths through the TC eye.
# this output looks nicer but kinda takes forever to run. Use xyplots_simple for a quicker function.
# Inputs: A dataframe holding shear corrected xy positions for different flight legs.
#         All flight legs provided will be added to the plot. Make sure to remove entries
#         before making the figure.
# Return: None. Side effect: make a flight path plot for the given passes.
def xyplots( df, limx=25, limy=25, title="Shear Relative Eye Flight Legs"):
    plt.figure(figsize = (5, 5))
    fs = 14
    scattersize = .2
    helper_fns_winter2023.change_font_sizes(fs,fs)
    plt.title(title)
    plt.ylabel("Upshear to Downshear Dist. (km)")
    plt.xlabel("Left to Right Shear Dist. (km)")
    plt.xlim([ -limx, limx])
    plt.ylim([ -limy, limy])

    # combine all x and y axes here
    xtotal, ytotal = [], []
    for index in range(len(df['xdistsshear'].values)):
        xtotal += df['xdistsshear'].values[index].tolist()
        ytotal += df['ydistsshear'].values[index].tolist()
    print("Total scatter points to add:" + str(len(xtotal)))

    dli, dri, uri, uli = 0, 0, 0, 0
    dlx, dly, drx, dry, urx, ury, ulx, uly = [], [], [], [], [], [], [], []
    for cloudi in range(len(xtotal)):
        xi, yi = xtotal[cloudi], ytotal[cloudi]

        # put x and y distances into correct quadrant lists
        # DL case
        if xi < 0.0 and yi > 0.0:
            dlx.append(xi)
            dly.append(yi)
        # DR case
        elif xi > 0.0 and yi > 0.0:
            drx.append(xi)
            dry.append(yi)
        # UR case
        elif xi > 0.0 and yi < 0.0:
            urx.append(xi)
            ury.append(yi)
        # UL case
        elif xi < 0.0 and yi < 0.0:
            ulx.append(xi)
            uly.append(yi)

    # do all the scatter plotting at once to save time!
    plt.scatter(dlx, dly, c='b', label="DL", s=scattersize * 20)
    plt.scatter(drx, dry, c='r', label="DR", s=scattersize * 20)
    plt.scatter(urx, ury, c='k', label="UR", s=scattersize * 20)
    plt.scatter(ulx, uly, c='y', label="UL", s=scattersize * 20)
    plt.legend(loc='upper right', framealpha=1, fontsize = fs * .6667)
    plt.show()


# the same as the function above, but it uses line plots instead of scatter plots to save plotting time
def xyplots_simple( df, limx=25, limy=25, title="Shear Relative Eye Flight Legs"):
    plt.figure(figsize = (5, 5))
    fs = 14
    scattersize = .2
    helper_fns_winter2023.change_font_sizes(fs,fs)

    plt.title(title)
    plt.ylabel("Upshear to Downshear Dist. (km)")
    plt.xlabel("Left to Right Shear Dist. (km)")
    plt.xlim([ -limx, limx])
    plt.ylim([ -limy, limy])


    for index in range(len(df['xdistsshear'])):
        plt.plot(df['xdistsshear'][index], df['ydistsshear'][index], c='k', linewidth=.5) # s=.25)

    plt.axvline(x=0.0, c='g', linewidth=2)
    plt.axhline(y=0.0, c='g', linewidth=2, label='Quadrants')
    plt.legend()