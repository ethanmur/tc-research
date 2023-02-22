# import modules
import numpy as np
import pandas as pd
from scipy.signal import find_peaks



# defintions of peak finder functions are saved outside of main code -> easy to switch in and out!
def find_peaks_simple( spd_avg):
    return find_peaks( x=spd_avg, height=40, prominence=10)[0]


# new method: find peaks using two slightly more inclusive methods, and take the union of the results!
# lower prominences but with tall heights! aysmmetric but strong tc
# or higher prominences with low heights

# this method doesn't work too well lol, even with the slightly higher prominences
# it's still finding too many false peaks
def find_peaks_union( spd_avg):
    vpeaks_lowprom = find_peaks( x=spd_avg, height=35, prominence=7)[0]
    vpeaks_lowh = find_peaks( x=spd_avg, height=25, prominence=10)[0]
    vpeaks = np.union1d(vpeaks_lowprom, vpeaks_lowh)

    return vpeaks


def find_peaks_pressure_mins( rawp, spd_avg, window):
    # find min pressure peaks, pmins
    p_avg = pd.Series( rawp).rolling(window=window, min_periods=1, center=True).mean()

    pmins_lowp = find_peaks( x= -p_avg, height=-960, prominence=10)[0]
    pmins_prominent = find_peaks( x= -p_avg, height=-990, prominence=20)[0]
    pmins = np.union1d( pmins_lowp, pmins_prominent)

    # print( 'pmins = ' + str( pmins))
    # print( 'length of spd_avg = ' + str( len( spd_avg)))

    # for now, just use the old find peaks method... update after testing if pmins works!

    # list holding all valid vpeaks. Two are generated for each pressure min
    vpeaks = []
    # find wind speed peaks to the left and right of each pressure minima
    for pmin in pmins:
        # print( 'pmin = ' + str(pmin))

        # speeds to the left of the pressure min
        spd_avg_left = spd_avg[ 0: pmin]
        # speeds to the right of the pressure min
        spd_avg_right = spd_avg[ pmin: len( spd_avg)]

        # find vmax values for this pressure min
        vpeaks_left = find_peaks( x=spd_avg_left, height=20, prominence=10)[0]
        vpeaks_right = find_peaks( x=spd_avg_right, height=20, prominence=10)[0]
        vpeaks_right = [ pmin + vpeak for vpeak in vpeaks_right ] # add by index to get back to correct values!

        # print( 'vpeaks left = ' + str( vpeaks_left))
        # print( 'vpeaks right = ' + str( vpeaks_right))

        # add the last left value and first right value to the vpeaks list!
        # closest values to pmin

        if len( vpeaks_left) != 0 and len( vpeaks_right) != 0:
            vpeaks.append( vpeaks_left[ -1])
            vpeaks.append( vpeaks_right[ 0])
        elif len( vpeaks_left) != 0:
            print( "No valid right peak found for this case!")
        elif len( vpeaks_right) != 0:
            print( "No valid left peak found for this case!")
        else:
            print( "No valid peaks found for this case!")

    # print( vpeaks)
    return vpeaks, pmins




# the same as the function above, but it limits the area around the low pressure center that
# can be searched! it finds the time at the center, and uses timelim as an upper and lower bound.
# default = 15 minutes on either side of the eye pass
def find_peaks_pmin_time_limit( rawp, spd_avg, time, window, timelim = 15, filter_inner_core=False):
    # find min pressure peaks, pmins
    # use two different, complementary methods
    p_avg = pd.Series( rawp).rolling(window=window, min_periods=1, center=True).mean()
    pmins_lowp = find_peaks( x= -p_avg, height=-960, prominence=10)[0]
    pmins_prominent = find_peaks( x= -p_avg, height=-990, prominence=20)[0]
    pmins = np.union1d( pmins_lowp, pmins_prominent)

    # test
    # print( 'pmins = ' + str( pmins))
    # print( 'length of spd_avg = ' + str( len( spd_avg)))

    # print( len( time))
    # print( time)

    # list holding all valid vpeaks. Two are generated for each pressure min
    vpeaks = []
    time_lims = []
    # find wind speed peaks to the left and right of each pressure minima
    for pmin in pmins:

        time_i = time[ pmin]

        # find indices for times to the left and right of center, closest to timelim!
        # / 60 to get to hours
        t_left_i = np.argmin( np.abs( time - ( time_i - timelim / 60) ) )
        t_right_i = np.argmin( np.abs( time - ( time_i + timelim / 60) ) )

        # 2/14/23 new code
        # if filter_inner_core has a non false value (would be a float wind speed),
        # check if there's a speed above the required threshold! if not, skip this
        # pmin
        if filter_inner_core:
            spd_avg_inner_core = spd_avg[ t_left_i : t_right_i]
            spd_max = np.nanmax( spd_avg_inner_core)

            # fail case: skip this pressure min
            if spd_max <= filter_inner_core:
                print("Curent pressure min fails this speed peak threshold!")
                continue

        # if no number was provided for filter_inner_core, or it passes the test
        # above, continue as expected!    
        time_lims.append( t_left_i)
        time_lims.append( t_right_i)
        # test
        # print( 'time limit = ' + str(timelim / 60))
        # print( 'pmin time = ' + str(time_i))
        # print( 'left time index = ' + str( time[t_left_i]))
        # print( 'right time index = ' + str( time[ t_right_i]))


        #####################
        # new, important code
        #####################

        # speeds to the left and right of the pressure min: within time limit!!
        spd_avg_left = spd_avg[ t_left_i : pmin]
        spd_avg_right = spd_avg[ pmin : t_right_i]


        # find vmax values for this pressure min
        # this code is more restrictive; first look for large eyewalls, then
        # look for smaller / asymmetric peaks

        vpeaks_left = find_peaks( x=spd_avg_left, height=25, prominence=10)[0]
        vpeaks_right = find_peaks( x=spd_avg_right, height=25, prominence=10)[0]


        # if any peaks are missing, use a more inclusive sorting algorithm
        if len( vpeaks_left) == 0 or len( vpeaks_right) == 0:
            vpeaks_left = find_peaks( x=spd_avg_left, height=15, prominence=5)[0]
            vpeaks_right = find_peaks( x=spd_avg_right, height=15, prominence=5)[0]

        # add really big peaks with basically no prominence
        if len( vpeaks_left) == 0 or len( vpeaks_right) == 0:
            vpeaks_left = find_peaks( x=spd_avg_left, height=40, prominence=1)[0]
            vpeaks_right = find_peaks( x=spd_avg_right, height=40, prominence=1)[0]


        # one last shot to add missing peaks! basically eliminate the prominence requirement
        if len( vpeaks_left) == 0 or len( vpeaks_right) == 0:
            vpeaks_left = find_peaks( x=spd_avg_left, height=12, prominence=2.5)[0]
            vpeaks_right = find_peaks( x=spd_avg_right, height=12, prominence=2.5)[0]


        # before adding to the list...
        # new code: need to correct left eyewall too! shift indices back to starting spot by t_left_i
        vpeaks_left = [ t_left_i + vpeak for vpeak in vpeaks_left ]
        vpeaks_right = [ pmin + vpeak for vpeak in vpeaks_right ] # add by index to get back to correct values!

        # testing
        # print( 'vpeaks left = ' + str( vpeaks_left))
        # print( 'vpeaks right = ' + str( vpeaks_right))

        # one last check, then add values to the eyewall list
        if len( vpeaks_left) != 0 and len( vpeaks_right) != 0:
            vpeaks.append( vpeaks_left[ -1])
            vpeaks.append( vpeaks_right[ 0])
        else:
            print( "Error in first if statement! there should alway be valid eyewalls")

    # print( vpeaks)
    return vpeaks, pmins, time_lims
