## last edited: 6/23/22
## Status: Automatically plots TDR and CRL cross sections with the same lat or lon coords.
##         This is working much better! some comparisons look perfect, others seem slightly off still (see 09/26 eye 1)

import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import datetime
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import make_plots

def choose_data( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here

    if tc.casefold() == 'fred':
        crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
        tdr_path = "/Users/etmu9498/research/data/tdr/fred/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_list = make_plots.load_crl(crl_path, print_files=False)

        crl_range = [
            (), (),
            (), (), (),
            (), (),
            (), ()]
        xlims = [
            ( -82.5, -68), ( -82.5, -68 ),
            ( -77, -72), ( -77, -72), (24, 18),
            ( -77, -73.5), (-77, -73.5 ),
            (-79, -75), (-79, -75) ]

        dates = [
            '08-11', '08-11',
            '08-12-am', '08-12-am', '08-12-am',
            '08-12-pm', '08-12-pm',
            '08-13', '08-13' ]
        eye_pass = ["1", "2", "1", "2", "3", "1", "2", "1", "2"]
        tc_name = 'Fred'

    elif tc.casefold() == 'grace':
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
        crl_range = []
        xlims = [
            (-85, -81), (-85, -81), (-85, -83), (24, 23), (-85, -83), (-85, -83),
            (-86.5, -82), (-88.5, -85.5), (-88, -86),
            (-90, -85), (-92, -87.5), (-92, -88.5), (-92, -88.5)]

        dates = [
            '08-27', '08-27', '08-27', '08-27', '08-27', '08-27', '08-27',
            '08-28', '08-28', '08-28', '08-29', '08-29', '08-29']
        eye_pass = [
            "1", "2", "3","4","5","6","7",
            "1", "2", "3", "1", "2", "3"]
        tc_name = 'Ida'

    elif tc.casefold() == 'ida-crl':
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

def choose_crl_date( date, crl_list):
    if date == '08-11':
        crl_data = crl_list[ 0]
    elif date == '08-12-am':
        crl_data = crl_list[ 1]
    elif date == '08-12-pm':
        crl_data = crl_list[ 2]
    elif date == '08-13':
        crl_data = crl_list[ 3]
    elif date == '08-16':
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
def plot( tc='sam', yaxis_zoom=False):
    warnings.filterwarnings("ignore")

    tcdata = choose_data( tc)
    if tcdata == 'selected TC name is not yet implemented':
        print( tcdata)
        return

    print( 'number of tdr files: ' + str( len( tcdata['tdr_list'] )))
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

        plt.figure(figsize=(18,12), facecolor='w')
        plt.subplot(311)
        plt.title( "TDR and CRL, Hurricane " + tcdata['tc_name'] + ", "
            + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])
        make_plots.plot_tdr( tcdata['tdr_path'], inbound_data, outbound_data, axis)
        plt.xlim( tcdata['xlims'] [counter][0], tcdata['xlims'] [counter][1])

        os.chdir( tcdata['crl_path'])
        data = xr.open_dataset( crl_data)

        if yaxis_zoom:
            # optional command to zoom in on TDR data only within CRL height range!
            plt.ylim( 0.0, np.max( - data.H))

        plt.subplot(312)
        make_plots.plot_power_ch1(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)
        plt.xlim( tcdata['xlims'][ counter][0], tcdata['xlims'][counter][1])


        # ## add a red line to power plot showing location when TDR scan was taken!
        # # get time from tdr header
        # float_time = float( inbound_data[9:11]) + float( inbound_data[11:13]) / 100
        #
        # for index in range( len( data.time.values)):
        #     time = data.time.values[index]
        #     # make sure the crl time and tdr scan time match to within .02 hours!
        #     if time % 1 <= float_time % 1 + .01 and time % 1 >= float_time % 1 - .01 and np.rint( time ) == np.rint( float_time):
        #         # print( float_time)
        #         # print( time)
        #         if axis == 'lat':
        #             plt.axvline( x= data.Lat[ index], c='r', linewidth=5)
        #         elif axis == 'lon':
        #             plt.axvline( x= data.Lon[ index], c='r', linewidth=5)
        #         break
        #

        plt.subplot(313)
        make_plots.plot_T(tcdata['crl_path'], crl_data, tcdata['crl_range'][counter][0], tcdata['crl_range'][counter][1], axis)
        plt.xlim( tcdata['xlims'][ counter][0], tcdata['xlims'][counter][1])


        if tc == 'ida-crl':
            os.chdir( "/Users/etmu9498/research/figures/tdr/ida/")
            plt.savefig( "tdr-crl.png" )
        elif yaxis_zoom:
            os.chdir( "/Users/etmu9498/research/figures/tdr/" + tcdata['tc_name'].casefold() + "/yzoom/")
            plt.savefig( "tdr-crl-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
        else:
            os.chdir( "/Users/etmu9498/research/figures/tdr/" + tcdata['tc_name'].casefold() + "/")
            plt.savefig( "tdr-crl-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
        print( "Image " + str( counter + 1) + " complete" )

        # this code makes a blank figure in between the real figures: it basically
        # just makes the output look prettier
        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(2) # figure height in inches

    warnings.filterwarnings("default")



def plot_tdr_only( tc='sam', yaxis_zoom=False):

    warnings.filterwarnings("ignore")

    tcdata = choose_data( tc)
    if tcdata == 'selected TC name is not yet implemented':
        print( tcdata)
        return

    print( 'number of tdr files: ' + str( len( tcdata['tdr_list'] )))
    for counter in range( len( tcdata[ 'dates'] )):


        # xlims are lats
        if tcdata[ 'xlims'] [ counter][0] > 0.0:
            axis = 'lat'
        # xlims are lons
        else:
            axis = 'lon'

        # load data
        os.chdir( tcdata['tdr_path'] )
        inbound_data = tcdata[ 'tdr_list'] [ counter*2]
        if len( tcdata['tdr_list']) > counter*2+1:
            outbound_data = tcdata['tdr_list'] [ counter*2 + 1]
        else:
            print( 'error!')

        plt.figure(figsize=(18,4), facecolor='w')
        plt.title( "TDR and CRL, Hurricane " + tcdata['tc_name'] + ", "
            + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter])
        make_plots.plot_tdr( tcdata['tdr_path'], inbound_data, outbound_data, axis)
        plt.xlim( tcdata['xlims'] [counter][0] , tcdata['xlims'] [counter][1] )

        os.chdir( tcdata['crl_path'])
        crl_data = choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )
        data = xr.open_dataset( crl_data)

        if yaxis_zoom:
            # optional command to zoom in on TDR data only within CRL height range!
            plt.ylim( 0.0, np.max( - data.H))
            os.chdir( "/Users/etmu9498/research/figures/tdr/" + tcdata['tc_name'].casefold() + "/tdr-only-yzoom/")
            plt.savefig( "tdr-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
        else:
            os.chdir( "/Users/etmu9498/research/figures/tdr/" + tcdata['tc_name'].casefold() + "/")
            plt.savefig( "tdr-" + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
        print( "Image " + str( counter + 1) + " complete" )

        # this code makes a blank figure in between the real figures: it basically
        # just makes the output look prettier
        f,ax = plt.subplots()
        f.set_visible(False)
        f.set_figheight(2) # figure height in inches

    warnings.filterwarnings("default")


# an old implementation of the plot() function, but without using a dictionary.
# I was checking to see if things ran faster this way, but they don't
# def plot_no_dict( tc='sam'):
#
#     warnings.filterwarnings("ignore")
#
#     tcdata = choose_data( tc)
#     crl_path = tcdata['crl_path']
#     tdr_path = tcdata['tdr_path']
#     crl_list = tcdata['crl_list']
#     tdr_list = tcdata['tdr_list']
#     crl_range = tcdata['crl_range']
#     xlims = tcdata['xlims']
#     dates = tcdata['dates']
#     eye_pass = tcdata['eye_pass']
#     tc_name = tcdata['tc_name']
#
#     print( 'number of tdr files: ' + str( len( tdr_list )))
#     for counter in range( len(crl_range )):
#
#         if xlims [ counter][0] > 0.0:
#             axis = 'lat'
#         else:
#             axis = 'lon'
#
#         # load data
#         os.chdir( tdr_path )
#         inbound_data = tdr_list [ counter*2]
#         outbound_data = tdr_list [ counter*2 + 1]
#         os.chdir( crl_path )
#
#         crl_data = choose_crl_date( dates[counter], crl_list )
#
#         plt.figure(figsize=(18,12), facecolor='w')
#         plt.subplot(311)
#         plt.title( "TDR and CRL, Hurricane " + tc_name + ", "
#             + dates[ counter] + ", Eye Pass " + eye_pass[ counter])
#         make_plots.plot_tdr(tdr_path, inbound_data, outbound_data, axis)
#         plt.xlim(xlims[counter][0], xlims[counter][1])
#
#         plt.subplot(312)
#         make_plots.plot_power_ch1( crl_path, crl_data, crl_range[counter][0], crl_range[counter][1], axis)
#         plt.xlim( xlims[ counter][0], xlims[counter][1])
#
#         ## add a red line to power plot showing location when TDR scan was taken!
#         # get time from tdr header
#         float_time = float( inbound_data[9:11]) + float( inbound_data[11:13]) / 100
#         data = xr.open_dataset( crl_data)
#         for index in range( len( data.time.values)):
#             time = data.time.values[index]
#             # make sure the crl time and tdr scan time match to within .02 hours!
#             if time % 1 <= float_time % 1 + .01 and time % 1 >= float_time % 1 - .01 and np.rint( time ) == np.rint( float_time):
#                 # print( float_time)
#                 # print( time)
#                 if axis == 'lat':
#                     plt.axvline( x= data.Lat[ index], c='r', linewidth=5)
#                 elif axis == 'lon':
#                     plt.axvline( x= data.Lon[ index], c='r', linewidth=5)
#                 break
#
#         plt.subplot(313)
#         make_plots.plot_T(crl_path, crl_data, crl_range[counter][0], crl_range[counter][1], axis)
#         plt.xlim( tcdata['xlims'][ counter][0], tcdata['xlims'][counter][1])
#
#
#         os.chdir( "/Users/etmu9498/research/figures/tdr/" + tc_name.casefold() + "/")
#         plt.savefig( "tdr-crl-" + tc_name.casefold() + "-" + str( counter+1) + ".png" )
#         print( "Image " + str( counter + 1) + " complete" )
#
#
#         f,ax = plt.subplots()
#         f.set_visible(False)
#         f.set_figheight(2) # figure height in inches
#
#     warnings.filterwarnings("default")
