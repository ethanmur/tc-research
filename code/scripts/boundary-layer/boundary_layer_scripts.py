import numpy as np
import os
import xarray as xr
import warnings

os.chdir("/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
os.chdir("/Users/etmu9498/research/code/scripts/in-situ-scripts")
import get_p3_heights







# create and return an xarray dataset with height corrected crl data for the whole dataset!
def height_correction_all_data( tcname, dataset, return_matrix='T_power'):

    metadata = tc_metadata.all_data(tc=tcname)

    crl_path = metadata['crl_path']
    crl_list = metadata['crl_list']
    crl_name = tc_metadata.choose_crl_date( metadata['dates'][dataset], crl_list)
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # get heights for flight level data!
    p3_heights = get_p3_heights.geth( tcname, dataset, all_data=True)

    if return_matrix == 'T_power':
        # rename variables
        T_2d = crl_data.T.where( crl_data.T.values < 50)
        power_2d = 10 * np.log10( crl_data.P_ch1 )
        power_2d = power_2d.where( power_2d.values > -30)

        # use a helper function to adjust for changing p-3 heights!
        newh, T_2d = get_p3_heights.interp_data2( T_2d, crl_data.H, p3_heights)
        print( 'temperature interpolated')
        newh, power_2d = get_p3_heights.interp_data2( power_2d, crl_data.H, p3_heights)
        print( 'power interpolated')

        return newh, T_2d, power_2d

    elif return_matrix == 'T_wv':
        # rename variables
        T_2d = crl_data.T.where( crl_data.T.values < 50)
        wv_2d = crl_data.WVMR.where( crl_data.WVMR.values != 0)
        wv_2d = wv_2d.where( wv_2d < 20)

        # use a helper function to adjust for changing p-3 heights!
        newh, T_2d = get_p3_heights.interp_data2( T_2d, crl_data.H, p3_heights)
        print( 'temperature interpolated')
        newh, wv_2d = get_p3_heights.interp_data2( wv_2d, crl_data.H, p3_heights)
        print( 'water vapor interpolated')

        return newh, T_2d, wv_2d



def find_pressure(tcname, dataset, slope = -10.257, height_correction=False):
    metadata = tc_metadata.all_data(tc=tcname)

    # load flight level and crl data
    flight_level_path = "/Users/etmu9498/research/data/in-situ-nc"
    flight_level_list = make_plots.load_flight_level( flight_level_path, print_files = False) # all in situ data names
    flight_level_name = tc_metadata.choose_in_situ_date( metadata['dates'][dataset], flight_level_list)
    os.chdir( flight_level_path)
    flight_data = xr.open_dataset( flight_level_name)

    crl_path = metadata['crl_path']
    crl_list = metadata['crl_list']
    crl_name = tc_metadata.choose_crl_date( metadata['dates'][dataset], crl_list)
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # make pressure matrix!!

    # special case to use height correction script automatically!
    if height_correction:
        print( 'correcting crl heights')
        h, temp, wv = height_correction_all_data( tcname, dataset, return_matrix='T_wv')
        height = np.nanmax( - h) * 1000

    # normal case: use regular crl data
    else:
        temp = crl_data.T
        power = crl_data.P_ch1


    # make a new matrix to hold pressure values; the same size as temp
    new_matrix=np.empty([ np.size( temp, 0), np.size( temp, 1)])

    # defining some constants
    height_count = np.shape( temp)[1]
    psurf = flight_data['PSURF.d'].values
    psurf = [ float( line) for line in psurf]

    warnings.filterwarnings("ignore")
    # do this for every timestep!

    print( 'creating pressure matrix')
    for crl_ind in range( np.shape( temp)[0]):

        # find the closest in situ timestep to current crl index
        # find closest crl distances to a 'padding' variable: save their indices and corresponding distances!
        insitu_idx = (np.abs( flight_data.time - crl_data.time[ crl_ind] )).argmin().values

        # get the surface pressure at the closest timestep
        psurfi = psurf[ insitu_idx]

        # make an equation relating pressure to height ( y = 0m here)
        yint = - ( slope * psurfi)

        if not height_correction:
            # get current p3 height
            height = float( flight_data['HT.d'][ insitu_idx].values)

        # find pressure at this height using the equation above!
        # height is defined in the if statement above
        ptop = ( height - yint) / slope

        # create a vertical pressure profile for this index
        pcolumn = np.linspace( ptop, psurfi, height_count)

        new_matrix[ crl_ind, :] = pcolumn


    warnings.filterwarnings("default")

    if height_correction:
        return h, temp, wv, new_matrix
    else:
        return new_matrix



def find_theta_thetav(tcname, dataset, slope = -10.257, height_correction=False):

    h, T_2d, wv_2d, pressure = find_pressure( tcname, dataset, height_correction=True)

    # calculate potential temperature!!! woo!!
    p0 = 1000
    r_cp = .286
    rl = 0
    theta = (T_2d + 273.15 ) * ( p0 / pressure)**( r_cp)
    thetav = theta * ( 1 + .61 * ( wv_2d / 1000) - rl  )

    return h, theta, thetav
