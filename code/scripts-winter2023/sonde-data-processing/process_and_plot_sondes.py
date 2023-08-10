# import...
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
from geopy import distance



# main function: loop through cases and create dropsonde plots! optionally, compare sonde results with CRL data.
def plot_all(tc='all', savefig=False):
    base_path = "/Users/etmu9498/research/data/dropsondes/ftp-nhc/"

    # option #1: do this for all sondes present
    if tc == 'all':
        os.chdir(base_path)
        file_names = []
        for (dirpath, dirnames, file) in os.walk( base_path):
            file_names.extend(dirnames)
            break

    # option #2: do this for one specified file
    elif type(tc) == type(""):
        file_names = [tc]	

    # option #3: do this for a list of files
    elif type(tc) == type([]):
        file_names = tc

    # do this for every folder in list (from 1 to all, depending on tc)
    for folderi, folder in enumerate( file_names):
        folder_path = base_path + "/" + folder + "/"

        print("Loading data for " + folder)

        # find the matching crl dataset name
        file_namesi = []
        for (dirpath, dirnames, file) in os.walk( folder_path):
            file_namesi.extend(file)
            break

        # save the original date here- helps figure out if rollover dates are present
        orig_date = ''
        # do this for every filename
        for filei, sounding_filename in enumerate( file_namesi):
            # get the original date:
            if filei == 0:
                orig_date = sounding_filename[1:9]

            # xarray loading doesn't work: use pandas scripts from bl final project instead!
            # example code found under ""
            col_names = [ 'IX', 't (s)', 'P (mb)', 'T (C)', 'RH (%)', 'Z (m)', 'WD', 'WS (m/s)', 'U (m/s)', 'V (m/s)', 'NS', 
                          'WZ (m/s)', 'ZW (m)', 'FP', 'FT', 'FH', 'FW', 'LAT (N)', 'LON (E)' ]

            # basically reading a text file here: skip header info (rows 1-21) and add automated column names (the same for every case??)
            # print(folder_path + sounding_filename)

            df = pd.read_fwf(folder_path + sounding_filename, skiprows=21, names=col_names)
            # give dataframe to a helper function for final processing
            process_and_plot( df, folder, sounding_filename, orig_date, filei, savefig)

# use this helper function to clean up sonde data
def process_and_plot( df, folder, sounding_filename, orig_date, count, savefig=False):
    # remove useless columns
    df = df.drop(columns=['IX', 'NS', 'FP', 'FT', 'FH', 'FW'])
    # replace missing values (currently 999's) with nans. Annoying placeholder format is fixed below.
    for col in df.columns:
        if col == 'Z (m)' or col == 'WD' or col == 'ZW (m)':
            df[col].replace(-999, np.nan,inplace=True)
        elif col == "P (mb)" or col == 'RH (%)' or col == 'WZ (m/s)':
            df[col].replace(-999.0, np.nan, inplace=True)
        if col == 'LAT (N)' or col == 'LON (E)':
            df[col].replace(-999.000, np.nan,inplace=True)
        else:
            df[col].replace(-999.00, np.nan,inplace=True)

    # remove nans for plotting purposes
    df = df.dropna(subset=('P (mb)', 'T (C)', 'RH (%)', 'Z (m)', 'WD', 'WS (m/s)'), how='any').reset_index(drop=True)

    # rename variables into more convenient forms and calculate a couple other variables
    T = df['T (C)'].values * units.degC
    z = df['Z (m)'].values * units.m
    p = df['P (mb)'].values * units.hPa
    wind_speed = df['WS (m/s)'].values * units.mps
    wind_dir = df['WD'].values * units.degrees
    RH = df['RH (%)'].values

    # use metpy for new calculations
    theta = mpcalc.potential_temperature(p, T)
    mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(p, T, RH / 100).to('g/kg')

    theta_v = mpcalc.virtual_potential_temperature(p, T, mixing_ratio)
    dewpoint = mpcalc.dewpoint_from_relative_humidity(T, RH * units.percent)
    theta_e = mpcalc.equivalent_potential_temperature(p, T, dewpoint)

    # get time info for nicer plotting
    date = sounding_filename[1:9]
    hours = sounding_filename[ 10:12]
    mins = sounding_filename[ 12:14]
    secs = sounding_filename[ 14:16]
    # if the current date doesn't match the original date, we have a rollover case! Aka a flight went from 23.99 hours to 24.01 hours.
    # in this case, add 24 to the time and use the original date
    if orig_date != date:
        time = 24 + float( hours) + float( mins) / 60 + float( secs) / 3600
        date = orig_date
    # normal case (no rollovers)
    else:
        time = float( hours) + float( mins) / 60 + float( secs) / 3600

    #print(orig_date)
    #print(date)
    #print()
    # getting distance from the TC center info here
    dist, crlmr, crltemp, crltheta, crlz, crltimestep = center_helper(date, time, p, z)

    # use the plotting helper function to look at the sonde data!
    plot_sonde_data( z, T, dewpoint, RH, theta, theta_v, theta_e, wind_speed, mixing_ratio, date, time, dist, count, crlmr = crlmr, crlt=crltemp, crltheta=crltheta, crlz=crlz, crltimestep=crltimestep)

    if savefig:
        path = "/Users/etmu9498/research/figures/dropsondes/"
        os.chdir( path)
        if not os.path.isdir( folder):
            os.makedirs( folder)
            print( 'New folder created: dropsondes/' + folder)
        savedir = path + folder
        os.chdir( savedir)
        plt.savefig( sounding_filename[:-3] + ".png", dpi=200, bbox_inches='tight')
        plt.close()
        print( "Plot " + sounding_filename + " saved" )



# make a simple 4 panel plot of sonde data, with optional crl data too!
# print and add some stats for convenience and comparison?
def plot_sonde_data( z, temp, dewpoint, RH, theta, theta_v, theta_e, winds, mr, date, time, dist, count, crlmr=False, crlt=False, crltheta=False, crlz=False, crltimestep=False):
    
    fig1, axs = plt.subplots(3, 3, gridspec_kw={'height_ratios': [1.5, 1.5, 1]}, figsize=(10, 10))
    

    # theta and theta e
    plt.subplot(331)    

    if len( z) > 0:
        plt.plot( theta, z, c='k', label='theta') 
        plt.plot( theta_v, z, c='y', label='theta_v')
        plt.plot( theta_e, z, c='r', label='theta_e')

        if type( crltheta) == type(np.array([])):
            plt.plot( crltheta, crlz, c='b', label='crl theta')
        plt.legend()
        plt.ylabel( "Sonde Height (m)")
        plt.xlabel( "Potential Temperatures (K)")
        # get limits for crl data plotting later!
        xlim = plt.gca().get_xlim()
        ylim = plt.gca().get_ylim()
        
        # relative humidity
        plt.subplot(332)
        plt.plot( RH, z, c='g') 
        plt.xlabel( "RH (%)")
        plt.ylabel("")
        plt.xlim( [60, 105])    
        
        # winds
        plt.subplot(333)
        plt.plot( winds, z, c='k') 
        plt.xlabel( "Wind Speed (m/s)")
        plt.ylabel("")

        
        # plot crl temp profiles!
        plt.subplot(334)
        plt.plot( dewpoint, z, c='b', label='Sonde T_d')    
        plt.plot( temp, z, c='r', label='Sonde T') 
        plt.ylabel( "Sonde Height (m)")
        plt.xlim([0, 35])

        if type(crlt) == type( np.array([])):
            plt.plot( crlt, crlz, c='k', label="CRL T")
        plt.legend()
        plt.xlabel( "Temperatures (K)")

        plt.text( 0, -1200, "Sonde Data: i = " + str(count) + ", date = " + date + ", time = " + str( time) + " hr, distance = " + str(dist) + ' km')
        if crltimestep:
            plt.text(0, -1600, 'time of crl sounding = ' + str(crltimestep) + ' hr')

        # plot mixing ratio
        plt.subplot(335)
        plt.plot( mr, z, c='Grey', label='Sonde MR')
        if type(crlmr) == type(np.array([])):
            plt.plot(crlmr, crlz, c='gold', label="CRL MR")
        plt.xlabel("WVMR (g/kg)")
        plt.ylabel("")
        plt.legend()

    axs[1, 2].axis('off')
    axs[2, 0].axis('off')
    axs[2, 1].axis('off')
    axs[2, 2].axis('off')


# pull center distances from crl data (much easier than using the original track datasets)
# date is in the following format: 20210929
# time is the decimal hours of the sonde launch
# along with returning the expected distances in km, return the temperature + height profile found at that level!
def center_helper(date, time, sonde_pressure, sonde_z):
    
    # go to the correct crl data folder
    crl_path = '/Users/etmu9498/research/data/crl-all-data-processed/' + date[0:4] + '/'
    os.chdir(crl_path)
    # find the matching crl dataset name
    file_names = []
    for (dirpath, dirnames, file) in os.walk( crl_path):
        file_names.extend(file)
        break
    
    crl_name = ''
    for name in file_names:
        # found the right name case!
        if date[4:9] == name[7:11]:
            crl_name = name
            break
    if crl_name == '':
        print("Error: no crl dataset found :/")
        return
    else:
        # load crl data
        crl_data = xr.open_dataset(crl_name)
    
    # get the closest time
    timeind = np.argmin(np.abs(crl_data.time.values - time))
    
    dist = crl_data.center_dist[timeind].values
    temp = crl_data.T[timeind, :].values
    mr = crl_data.WVMR[timeind, :].values
    z = crl_data.height.values

    # calculate potential temperature
    theta = find_theta(crl_data, sonde_pressure, sonde_z, timeind)

    return dist, mr, temp, theta, z, crl_data.time[timeind].values

# using the pressure field from the dropsonde and the crl temperature field, calculate crl theta!
def find_theta(crl_data, sonde_pressure, sonde_z, crli):
    # locally define crl height and temp
    crltemp = crl_data.T[crli, :].values
    crlz = crl_data.height.values
    
    if len(crltemp) == 0 or len(sonde_z) == 0:
        return( np.array([]))

    p = []
    # get the closest pressure values for every crl height level!
    for hi, hval in enumerate(crlz):
        # find the index of the closest sonde height
        sondei = np.argmin(np.abs( np.array(sonde_z) - hval))

        # if the difference is too large, just append a nan
        if np.nanmin(np.abs( np.array(sonde_z) - hval)) > 8:
            p.append(np.nan)
        # the heights are close to one another! no errors
        else:
            p.append( sonde_pressure[sondei].magnitude)
    p = np.array(p) * units.hPa
    crltemp = crltemp * units.degC
    theta = mpcalc.potential_temperature(p, crltemp)

    return theta.magnitude