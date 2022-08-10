# the two functions below automate the colorbar and line plot creation processes
# and automatically save figures for all, or a specified, tc to a folder in the
# figures tab.
# edited 8/9/22

import os
import matplotlib.pyplot as plt
import warnings

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
os.chdir(  "/Users/etmu9498/research/code/scripts/in-situ-scripts")
import in_situ_colorbar_lines

def flight_level_colorbar_auto( tc='all', variable_list=['rr']):

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    for tcname in tcname_list:
        warnings.filterwarnings("ignore")

        # maybe switch to all data here!
        tcdata = tc_metadata.choose_data_cloud_tops_good_data( tcname)
        if tcdata == 'selected TC name is not yet implemented':
            print( tcdata)
            return

        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['xlims'] ))+ '\n')

        for counter in range( len( tcdata[ 'dates'] )):

            if tcdata[ 'xlims'] [ counter][0] > 0.0:
                axis = 'lat'
            else:
                axis = 'lon'

            # load data
            os.chdir( tcdata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

            os.chdir( tcdata['in_situ_path'] )
            in_situ_data = tc_metadata.choose_in_situ_date( tcdata['dates'][counter], tcdata['in_situ_list'] )

            # not sure why this is commented out
            # make a crl plot with variable colorbars for all tc cases
            '''
            if len ( tcdata['crl_range'][counter]) == 2:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][1]
            elif len( tcdata['crl_range'][counter]) == 4:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][3]
            else:
                print( 'Error in number of indices! Update them in tc_metadata.py')
            '''

            xlims = [ tcdata['xlims'][counter][0], tcdata['xlims'][counter][1] ]

            title = "CRL and In Situ Data, TC " + tcdata['tc_name'] + ", " + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter]
            in_situ_colorbar_lines.flight_level_colorbar( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_data, tcdata['crl_range'][counter], axis, variable_list, xlims, title)

            # case to plot and save lat or lon values vs time!
            # helpful for determining if i1 and i2 are accurate
            # comment out some of the matching code above to get this to work
            # title = "In Situ Lat Lon Peaks, TC " + tcdata['tc_name'] + ", " + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter]
            # plot_lat_lon( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_data, i1, i2, axis, variable_list, xlims, title)
            # os.chdir( "/Users/etmu9498/research/figures/lat-lons/")

            # save the figure!
            os.chdir( "/Users/etmu9498/research/figures/in-situ-colorbars/")
            fig_name = ''
            for in_situ_var_name in variable_list:
                fig_name += in_situ_var_name
                fig_name += '-'

            plt.savefig( fig_name + tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png", bbox_inches='tight')
            print( "Image " + str( counter + 1) + " complete" )



def flight_level_lines_auto( tc='all', variable_list=['uwz', 'ws', 'rr']):

    # put tcname into a list to make the for loop work correctly
    if tc == 'all':
        tcname_list = ['grace', 'henri', 'ida', 'sam']
    else:
        tcname_list = [ tc]

    for tcname in tcname_list:
        warnings.filterwarnings("ignore")

        tcdata = tc_metadata.choose_data_cloud_tops_good_data( tcname)
        if tcdata == 'selected TC name is not yet implemented':
            print( tcdata)
            return

        print( "\nTC " + tcdata['tc_name'])
        print( 'Number of crl files: ' + str( len( tcdata['xlims'] ))+ '\n')

        for counter in range( len( tcdata[ 'dates'] )):

            if tcdata[ 'xlims'] [ counter][0] > 0.0:
                axis = 'lat'
            else:
                axis = 'lon'

            # load data
            os.chdir( tcdata['crl_path'] )
            crl_data = tc_metadata.choose_crl_date( tcdata['dates'][counter], tcdata['crl_list'] )

            os.chdir( tcdata['in_situ_path'] )
            in_situ_data = tc_metadata.choose_in_situ_date( tcdata['dates'][counter], tcdata['in_situ_list'] )

            '''
            # make a crl plot with variable colorbars for all tc cases
            if len ( tcdata['crl_range'][counter]) == 2:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][1]
            elif len( tcdata['crl_range'][counter]) == 4:
                i1 = tcdata[ 'crl_range'][counter][0]
                i2 = tcdata[ 'crl_range'][counter][3]
            else:
                print( 'Error in number of indices! Update them in tc_metadata.py')
            '''

            xlims = [ tcdata['xlims'][counter][0], tcdata['xlims'][counter][1] ]
            title = "CRL and In Situ Data, TC " + tcdata['tc_name'] + ", " + tcdata['dates'][ counter] + ", Eye Pass " + tcdata['eye_pass'][ counter]
            in_situ_colorbar_lines.flight_level_lines( tcdata['crl_path'], crl_data, tcdata['in_situ_path'], in_situ_data, tcdata['crl_range'][counter], axis, variable_list, xlims, title)

            # save the figure!
            os.chdir( "/Users/etmu9498/research/figures/in-situ-lines/")
            plt.savefig( tcdata['tc_name'].casefold() + "-" + str( counter+1) + ".png" )
            print( "Image " + str( counter + 1) + " complete" )
