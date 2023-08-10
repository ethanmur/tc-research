# import...
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import metpy.calc as mpcalc
from metpy.units import units
from geopy import distance
import scipy

# load some of the helper functions from the other dropsonde script
os.chdir("/Users/etmu9498/research/code/scripts-winter2023/sonde-data-processing")
import process_and_plot_sondes

def find_errors(tc='all'):
    base_path = "/Users/etmu9498/research/data/dropsondes/ftp-nhc/"
    file_names = get_folders(tc, base_path)

    filecount = 0

    # make a pandas dataframe to hold crl vs sonde values for each sonde! If there are 700 sondes, there will be 700 df rows
    # Will have the following structure:
    # year   |   date   | sonde #  | heights  | radial dist | crlt   |  sondet
    # -----------------------------------------
    #  ...   |    ...   |   ....   |    ...   |  ...        | ...    |  ...
    df_peaks = pd.DataFrame( )
    # save values in temporary lists here: will be added to df_peaks once filled
    df_yearlist = []
    df_datelist = []
    df_sondecountlist = []
    df_timelist = []
    df_heightslist = []
    df_distlist = []
    df_crlt = []
    df_sondet = []
    df_crlwv = []
    df_sondewv = []


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

            # getting distance from the TC center + thermodynamic comparison info here 
            # define locally to match older code, taken from "process_and_plot_sondes.center_helper()"
            sonde_pressure, sonde_z = p, z
            
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
            
            # get the closest crl time for this drop
            timeind = np.argmin(np.abs(crl_data.time.values - time))
            
            # get info from the closest crl profile
            radial_dist = crl_data.center_dist[timeind].values
            crlt = crl_data.T[timeind, :].values
            crlwv = crl_data.WVMR[timeind, :].values
            crlz = crl_data.height.values
            # calculate potential temperature using a helper function
            crltheta = process_and_plot_sondes.find_theta(crl_data, sonde_pressure, sonde_z, timeind)

            # interpolate dropsonde temp, wv, and theta down to crl grid spacing
            if len(crlt) == 0 or len(z) == 0:
                sondet = np.array([])
                sondewv = np.array([])
                sondetheta = np.array([])
            else:

                # find the first (highest) height that has non-nan values!
                heightinds = np.where( ~ np.isnan( crlt))[0]
                if len(heightinds) == 0:
                    empty = np.array([])
                    crlz, crlt, crlwv, sondet, sondewv = empty, empty, empty, empty, empty
                else:
                    crlz = crlz[ heightinds[0] :]
                    crlt = crlt[ heightinds[0] :]
                    crlwv = crlwv[ heightinds[0] :]

                    # the numpy version doesn't interpret above bounds correctly :/
                    # flip all arrays because np.interp() needs increasing x axes
                    sondet = np.interp( np.flip(crlz), np.flip(z.magnitude), np.flip(T.magnitude))
                    sondewv = np.interp( np.flip(crlz), np.flip(z.magnitude), np.flip(mixing_ratio.magnitude))
                    sondet = np.flip(sondet)
                    sondewv = np.flip(sondewv)
                    # sondetheta = np.interp( crlz, z.magnitude, theta.magnitude)

                    # trying to use scipy for the interpolation... having issues with nan inputs :/ using np method instead!
                    '''
                    sondetf = scipy.interpolate.interp1d( crlz, crlt)
                    sondewvf = scipy.interpolate.interp1d( crlz, crlwv)
                    sondethetaf = scipy.interpolate.interp1d( crlz, crltheta)
                    sondet = sondetf( T.magnitude)
                    sondewv = sondewvf( mixing_ratio.magnitude)
                    sondetheta = sondethetaf( theta.magnitude)
                    '''

                    '''
                    np.set_printoptions(threshold=np.inf)                
                    print(crlz)
                    print(z.magnitude)
                    print(crlt)
                    print(T.magnitude)
                    print(sondet)
                    np.set_printoptions(threshold=1000)
                    '''

            # add sonde and crl info to the dataframe!
            df_yearlist.append(sounding_filename[1:5])
            df_datelist.append(date[4:])
            df_sondecountlist.append(filei)
            df_timelist.append(crl_data.time[timeind].values)
            df_heightslist.append(crlz)
            df_distlist.append(radial_dist)
            df_crlt.append(crlt)
            df_sondet.append(sondet)
            df_crlwv.append(crlwv)
            df_sondewv.append(sondewv)

            # iterate number of sondes... nice for testing
            filecount += 1
        
    print( "Number of Sondes Used: " + str(filecount))

    df_peaks['year'] = df_yearlist
    df_peaks['date'] = df_datelist
    df_peaks['sonde number'] = df_sondecountlist
    df_peaks['time (UTC)'] = df_timelist
    df_peaks['height (m)'] = df_heightslist
    df_peaks['radial distance (km)'] = df_distlist
    df_peaks['crl t (C)'] = df_crlt
    df_peaks['sonde t (C)'] = df_sondet
    df_peaks['crl wv (g/kg)'] = df_crlwv
    df_peaks['sonde wv (g/kg)'] = df_sondewv

    return df_peaks


def get_folders(tc, base_path):
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
    return file_names


# make simple plots comparing crl and sonde t and wv errors
# the first two subplots just plot each crl sonde comparison as a scatter point,
# while the second two create a colormap of point probabilities.
def simple_plot(df=False, height_limit=[], bincount=50):

    if len(height_limit) == 2:
        allcrlt, allcrlwv, allsondet, allsondewv, tanom, wvanom, theights, wvheights = plot_helper(df, height_limit)
    else:
        allcrlt, allcrlwv, allsondet, allsondewv, tanom, wvanom = plot_helper(df, height_limit)
    
    plt.figure(figsize=(10, 10))

    plt.subplot(221)
    plt.scatter(allcrlt, allsondet, facecolors='none', edgecolors='g', s=10, marker='o')
    plt.plot( np.arange(-10, 35), np.arange(-10, 35), c='k', linewidth=.75)
    plt.xlabel("CRL T (C)")
    plt.ylabel("Sonde T (C)")
    plt.title("Temperature Comparisons")
    plt.xlim([-10, 35])
    plt.ylim([-10, 35])

    plt.subplot(222)
    plt.scatter(allcrlwv, allsondewv, facecolors='none', edgecolors='r', s=10, marker='o')
    plt.plot( np.arange(-10, 35), np.arange(-10, 35), c='k', linewidth=.75)
    plt.xlabel("CRL WVMR (g/kg)")
    plt.ylabel("Sonde WVMR (g/kg)")
    plt.title("WVMR Comparisons")
    plt.xlim([-10, 30])
    plt.ylim([-10, 30])

    map = plt.cm.get_cmap( "RdYlBu").reversed()
    xmin, xmax, ymin, ymax = -10, 35, -10, 35
    plt.subplot(223)
    plt.hist2d(allcrlt, allsondet, cmap=map, bins=bincount, range=[[xmin, xmax], [ymin, ymax]], norm=matplotlib.colors.LogNorm()) # int(len(allcrlt)/1000))
    plt.plot( np.arange(xmin, xmax), np.arange(xmin, xmax), c='k', linewidth=.75)
    plt.xlabel("CRL T (C)")
    plt.ylabel("Sonde T (C)")
    plt.title("Temperature Comparisons")
    plt.colorbar(label="")
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])

    plt.subplot(224)
    xmin, xmax, ymin, ymax = -10, 30, -10, 30
    plt.hist2d(allcrlwv, allsondewv, cmap=map, bins=bincount, range=[[xmin, xmax], [ymin, ymax]], norm=matplotlib.colors.LogNorm()) # int(len(allcrlwv)/1000))
    plt.plot( np.arange(xmin, xmax), np.arange(xmin, xmax), c='k', linewidth=.75)
    plt.xlabel("CRL WVMR (g/kg)")
    plt.ylabel("Sonde WVMR (g/kg)")
    plt.title("WVMR Comparisons")
    plt.colorbar(label="Bin Count")
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])

# plot crl - sonde anomaly histograms for t and wv here
# also plot general crl and sonde distributions for comparison
def anom_plot(df, height_limit=[], xmin=-10, xmax=35, bincount=50, bw=5):
    if len(height_limit) == 2:
        allcrlt, allcrlwv, allsondet, allsondewv, tanom, wvanom, theights, wvheights = plot_helper(df, height_limit)
    else:
        allcrlt, allcrlwv, allsondet, allsondewv, tanom, wvanom = plot_helper(df, height_limit)

    print("Temp Anomaly Mean = " + str(np.mean(tanom)))
    print("Temp Anomaly std = " + str(np.std(tanom)))
    print("WV Anomaly Mean = " + str(np.mean(wvanom)))
    print("WV Anomaly std = " + str(np.std(wvanom)))

    xinc = np.arange(xmin,xmax,bw)
    colors = ['b', 'c', 'skyblue', 'r', 'g', 'lightgreen']

    # min4hist=np.round(np.nanmin( data_list[ i]),1)-binsize[ i]
    # max4hist=np.round(np.nanmax( data_list[ i]),1)+binsize[ i]
    # nbins=int((max4hist-min4hist)/binsize[ i])
    # plt.hist( data_list[ i],nbins,edgecolor='black', color=colors[ i])
    
    # special height limit case: print info confirming the tests worked
    if len(height_limit) == 2:
        print("maxs")
        print(np.nanmax(np.array(theights)))
        print(np.nanmax(np.array(wvheights)))
        print('mins')
        print(np.nanmin(np.array(theights)))
        print(np.nanmin(np.array(wvheights)))

    plt.figure(figsize=(12, 8))

    plt.subplot(231)
    plt.title("CRL T (C)")
    plt.ylabel("Count")
    hx=np.histogram( allcrlt, xinc)
    plt.bar(hx[1][:-1], hx[0], edgecolor='k', color=colors[0], width=bw, linewidth=1)
    plt.xlim([xmin, xmax])

    plt.subplot(232)
    hx=np.histogram( allsondet, xinc)
    plt.bar(hx[1][:-1], hx[0], edgecolor='k', color=colors[1], width=bw, linewidth=1)
    plt.title("Sonde T (C)")
    plt.xlim([xmin, xmax])

    plt.subplot(233)
    plt.title("T Anomaly (CRL - Sonde)")
    hx=np.histogram( tanom, xinc)
    plt.bar(hx[1][:-1], hx[0], edgecolor='k', color=colors[2], width=bw, linewidth=1)
    plt.xlim([xmin, xmax])
    plt.axvline(x=0, c='k', linewidth=.75)
    
    plt.subplot(234)
    plt.title("CRL WVMR (g/kg)")
    hx=np.histogram( allcrlwv, xinc)
    plt.bar(hx[1][:-1], hx[0], edgecolor='k', color=colors[3], width=bw, linewidth=1)
    plt.xlim([xmin, xmax])
    plt.ylabel("Count")

    plt.subplot(235)
    plt.title("Sonde WVMR (g/kg)")
    hx=np.histogram( allsondewv, xinc)
    plt.bar(hx[1][:-1], hx[0], edgecolor='k', color=colors[4], width=bw, linewidth=1)
    plt.xlim([xmin, xmax])

    plt.subplot(236)
    plt.title("WVMR Anomaly (CRL - Sonde)")
    hx=np.histogram( wvanom, xinc)
    plt.bar(hx[1][:-1], hx[0], edgecolor='k', color=colors[5], width=bw, linewidth=1)
    plt.xlim([xmin, xmax])
    plt.axvline(x=0, c='k', linewidth=.75)
    

# helper function used to simplify the two plotting scripts above.
# return all valid crl + sonde points here. 
# also return t and wv anomalies (calculated from each point, not subtracting overall distributions later)
def plot_helper(df, height_limit):
    # no crl sonde dataframe provided case
    if type(df) != type(pd.DataFrame()):
        df = find_errors(tc='all')

    # save valid crl and sonde combos here
    allcrlt, allcrlwv = [], []
    allsondet, allsondewv = [], []

    # save anomalies (crl - sonde values) here
    tanom, wvanom = [], []
    theights, wvheights = [], []

    # loop through each row aka unique sonde
    for rowi in range(len(df['year'])):        

        crlt = df['crl t (C)'][rowi]
        sondet = df['sonde t (C)'][rowi]        
        crlwv = df['crl wv (g/kg)'][rowi]
        sondewv = df['sonde wv (g/kg)'][rowi]        

        # find valid data point inds
        nonan_inds = np.intersect1d( np.where( ~np.isnan(crlt))[0], np.where( ~np.isnan(sondet))[0])
        # check if we need to sort by height here!
        # first index is the lower bound, second is the upper
        if len(height_limit) == 2:
            nonan_inds2 = np.intersect1d( np.where( df['height (m)'][rowi] >= height_limit[0])[0], np.where( df['height (m)'][rowi] < height_limit[1])[0])
            nonan_inds = np.intersect1d( nonan_inds, nonan_inds2)
            theights += df['height (m)'][rowi][nonan_inds].tolist()

        allcrlt += crlt[nonan_inds].tolist()
        allsondet += sondet[nonan_inds].tolist()

        # calculate and save temperature anomalies
        tanom += (crlt[nonan_inds] - sondet[nonan_inds]).tolist()

        # repeat for water vapor measurements here
        nonan_inds = np.intersect1d( np.where( ~np.isnan(crlwv))[0], np.where( ~np.isnan(sondewv))[0])
        if len(height_limit) == 2:
            nonan_inds2 = np.intersect1d( np.where( df['height (m)'][rowi] >= height_limit[0])[0], np.where( df['height (m)'][rowi] < height_limit[1])[0])
            nonan_inds = np.intersect1d( nonan_inds, nonan_inds2)            
            wvheights += df['height (m)'][rowi][nonan_inds].tolist()

        allcrlwv += crlwv[nonan_inds].tolist()
        allsondewv += sondewv[nonan_inds].tolist()

        # calculate and save temperature anomalies
        wvanom += (crlwv[nonan_inds] - sondewv[nonan_inds]).tolist()

    # make the plots!
    print("Number of valid T comparisons: " + str(len(allcrlt)))
    print("Number of valid WV comparisons: " + str(len(allsondewv)))

    if len(height_limit) == 2:
        return allcrlt, allcrlwv, allsondet, allsondewv, tanom, wvanom, theights, wvheights
    else:
        return allcrlt, allcrlwv, allsondet, allsondewv, tanom, wvanom