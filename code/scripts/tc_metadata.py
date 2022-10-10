# a file containing x axes, dates, and other metadata for plotting tc data automatically!
# other python scripts call these functions to load in data as dictionaries
import os
os.chdir( "/Users/etmu9498/research/code/scripts/")
import make_plots

# main function holding metadata for each TC!
# when called, this function returns a dictionary of values for each respective TC.
# cases for complete datasets only are included.
def all_data( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here
    # paths to data
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    in_situ_path = "/Users/etmu9498/research/data/in-situ"

    new_tdr_path = "/Users/etmu9498/research/data/tdr-new"
    new_crl_path = "/Users/etmu9498/research/data/crl-new"
    new_flight_data_path = '/Users/etmu9498/research/data/in-situ-new'

    updated_matrices_crl_path = "/Users/etmu9498/research/data/crl-new-matrices"

    # load a list of available data to these variables
    crl_list = make_plots.load_crl(crl_path, print_files=False)
    in_situ_list = make_plots.load_flight_level(in_situ_path, print_files=False)

    if tc.casefold() == 'fred':
        tdr_path = "/Users/etmu9498/research/data/tdr/fred/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        xlims = [
            ( -77, -72), ( -77, -72), (24, 18),
            ( -77, -73.5), (-77, -73.5 ),
            (-79, -75), (-79, -75) ]
        dates = [
            '08-12-am', '08-12-am', '08-12-am',
            '08-12-pm', '08-12-pm',
            '08-13', '08-13' ]
        eye_pass = [ "1", "2", "3", "1", "2", "1", "2"]
        tc_name = 'Fred'

    elif tc.casefold() == 'grace':
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        crl_range = [ ( 0, 1800), (2250, 4100), (4800, 6550), # 8/17
            (0, 2100), (2800, 3570, 3655, 4400), (5000, 6200), # 8/18
            # eye 2 had a spiral, 4 indices get rid of it
            (0, 1600), (2100, 3300), (4000, 5700) ] # 8/19 # eye 2 (2100, 3900) removed... add it again!

        xlims = [ (20, 16), (-78.5, -74), (20, 17), # 8/17
            ( 22, 18), (-86.5, -82.5 ), (21, 18 ), # 8/18
            (-92, -89), (-92, -90), (-93, -89) ] # 8/19
        xtype = ['lat', 'lon', 'lat', 'lat', 'lon', 'lat', 'lon', 'lon', 'lon']

        # original eyewall distances when using tdr axis
        eyewall_dists = [ ( -12.5, 17.5), (-35, 12.5  ), (-10, 35),
            (-62.5, -10), (10, 85), ( -55, 25) ]
        # new eyewall distances for height corrected data with an x axis based on in situ data
        # new, precise limits
        in_situ_eyewall_dists = [ (-65, 72.5), (11, 57), (-12.5, 53), # old pass 2 and 3 eyewalls: (-75, 57.5), (-48, 53),
                (-11.5, 17), (-22.0, 7.5), (-12.5, 26),
                (-32.5, 11), ( 10, 70), (-55.5, 29.5)]

        # old limits
        # [ (-12.5, 15), (-25, 10), (-20, 25),
        #     (-35, 12.5), (10, 70), (-55, 30)]


        dates = [
            '08-17', '08-17', '08-17',
            '08-18', '08-18', '08-18',
            '08-19', '08-19', '08-19' ]
        eye_pass = [ "1", "2", "3", "1", "2", "3", "1", "2", "3"]

        # parameters used for getting an accurate distance axis in new crl data
        stretch = [ 1, 1, 1, 1, 1, 1]
        shift = [ 20, 0, 0, 0, 0, 10]

        tc_name = 'Grace'

        shrd = [ 107, 107, 107, 105, 105, 105, 81, 81, 81] # 0 utc the following day used for most cases
        shtd = [ 142, 142, 142, 135, 135, 135, 169, 169, 169]

        shear_quads = [ ('UL', "DR"), ('DL', 'UR'), ('DL', 'UR'),
                    ('UR', 'DL'), ('DL', 'UR'), ('DR', 'UL'),
                    ('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL')] # determined with shear_angles.py scripts! saved here for convenience

        shear_dir = shtd
        shear_mag = shrd

        crl_shear_inds = []
        tc_center_lat = []
        tc_center_lon = []

        # tc intensity (kts) at time of flight, taken from noaa summary docs
        intensity = [45, 45, 45, 70, 70, 70, 55, 55, 55 ]
        intensity_cat = find_cat( intensity)

    elif tc.casefold() == 'henri':
        tdr_path = "/Users/etmu9498/research/data/tdr/henri/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_range = [
            # add at least eye 1 and 3 from 8/20 dataset! Even though the clouds aren't
            # too tall, there is a distinct warm core in all 3 eye passes. Pass 2
            # is missing a bit of crl data though :(
            ( 400, 1800), ( 4970, 6050), # 8/20
            (0, 1565, 1650, 2400), (3100, 5100), (6700, 8100) ] # 8/21
            # eye 1 had a spiral, 4 indices get rid of it

        xlims = [
            ( 33.5, 29.5), ( -75, -72),
            ( 39, 33.5), (-73, -67.5), (-72, -68.5) ]
        xtype = [ 'lat', 'lon', 'lat', 'lon', 'lon']

        eyewall_dists = [(-22, 57), (-30, 75),
            (-12.5, 35), (-10, 30), (-20, 30)]

        # a lot of these values need to be improved :(
        in_situ_eyewall_dists = [( 5, 45), (5, 36), # fix only positive dists!
            (-11, 21), (-17.5, 25), (-21, 36)]

            # [(-25, 55), (-30, 77.5),
            # (-10, 20), (-17.5, 25), (-20, 35)]


        dates = [ '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        stretch = [ 1, 1, 1, 1, 1] # .66667
        shift = [ -50, 0, 0, 0, -5] # -12.5

        eye_pass = [ "1", "3", "1", "2", "3"]
        tc_name = 'Henri'

        shrd = [ 128, 128, 60, 60, 60]
        shtd = [ 179, 179, 240, 240, 240]
        shear_dir = shtd
        shear_mag = shrd
        shear_quads = [ ('DR', 'UL'), ('UL', 'DR'),
                ('DL', 'UR'), ('DR', 'UL'), ('UL', 'DR')]


        crl_shear_inds = []
        tc_center_lat = []
        tc_center_lon = []

        intensity = [60, 60, 65, 65, 65 ]
        intensity_cat = find_cat( intensity)

    elif tc.casefold() == 'ida':
        tdr_path = "/Users/etmu9498/research/data/tdr/ida/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        # list of touples representing the min and max indices used to clip crl data
        crl_range = [ ( 1000, 2600), (4800, 6400), (10200, 12000), (1500, 4700)]
        xlims = [ (-85, -81), (-85, -82), (-86, -82), (-92, -87)]
        xtype = [ 'lon', 'lon', 'lon', 'lon']

        # maybe add data from early 8/28? Looks like the P-3 was sampling something
        # else, but worth a closer look... maybe 1 good eye dataset here!

        eyewall_dists = [(-15, 65), (-35, 15), (-30, 50), (-12.5, 10)]
        in_situ_eyewall_dists = [(-25, 60), (-26, 34), (-15, 75), (-7, 20)]

        # older in situ eyewalls: new values are more accurate
        # [(-25, 65), (-25, 35), (-15, 75), (-25, 20)]

        dates = [ '08-27', '08-27', '08-27', '08-29']
        stretch = [ 1, 1, 1, 1]
        shift = [ 0, 0, 0, 90]

        eye_pass = [ "1", "2", "7", "2"]
        tc_name = 'Ida'

        shrd = [ 115, 115, 115, 112] # 18 utc used for last case
        shtd = [ 61, 61, 61, 127]

        crl_shear_inds = []

        shear_dir = shtd
        shear_mag = shrd
        shear_quads = [ ('UL', 'DR'), ('DL', 'UR'), ('DL', 'UR'), ('DL', 'UR')]

        tc_center_lat = []
        tc_center_lon = []

        intensity = [70, 70, 70, 130 ]
        intensity_cat = find_cat( intensity)

    elif tc.casefold() == 'sam':
        tdr_path = "/Users/etmu9498/research/data/tdr/sam/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        # list of touples representing the min and max indices used to clip crl data
        crl_range = [
            (0, 1900), (2300, 3800), (4400, 5600),
            (2600, 3900), (4800, 6500), # took eye 1 (0, 2000) out... only 1/2 the crl data is there :(
            (0, 2000), (2100, 3800)]
        xlims = [
            (-52, -49), (-51.5, -49.75), (-51.75, -50),
            (18, 15), (-54, -52), # (-55, -51),
            (-59.5, -56.5), (-59.5, -56.5) ]
        xtype = [ 'lon', 'lon', 'lon', 'lat', 'lon', 'lon', 'lon']

        eyewall_dists = [( -5, 15), (-10, 15), (-30, 12.5),
            (-5, 30), (-20, 10),
            (-25, 20), (-15, 35) ]

        # case 9/26 eye 3 is hard here: wind speed max is around -40, but cloud towers continue to around -47 km...
        in_situ_eyewall_dists = [ (-3, 31), (-2, 28), (-40, 5),
                # case 2 is hard... moved eyewall back to match plume
                (-2, 37), (-22, 12.5),
                (-33, 33), (-21, 34)]

        # older in situ eyewalls: new values are more accurate
        # [ (-5, 30), (-5, 30), (-40, 5),
        #     (-5, 35), (-20, 12.5), (-35, 30), (-20, 35)]

        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-29', '09-29']
        stretch = [ 1, 1, 1, 1, 1, 1, 1]
        shift = [ 65, 0, 0, 0, 0, 0, 0]

        eye_pass = [ "1", "2", "3", "2", "3", "1", "2"]
        tc_name = 'Sam'
        # shear direction taken at 0 UTC the following day for each case
        # Ex: sam 9/26 passes ran from 22 to 25 UTC, closest shear values at 0
        # UTC on 9/27

        # possible value headers:
        # SHDC: : Same as SHRD but with vortex removed and averaged from 0-500 km relative
        # to 850 hPa vortex center
        # SDDC: heading of SHDC. Westerly shear == 90 degrees
        shdc = [78, 78, 78, 54, 54, 56, 56 ]
        sddc = [48, 48, 48, 12, 12, 25, 25 ]

        # SHRD: 850-200 hPa shear magnitude (kt *10) vs time (200-800 km)
        # SHTD: heading of SHRD
        shrd = [ 90, 90, 90, 85, 85, 103, 103 ]
        shtd = [ 64, 64, 64, 44, 44, 57, 57]

        # SHRS: 850-500 hPa shear magnitude (kt *10) vs time
        # SHTS: heading of SHRS
        shrs = [ 33, 33, 33, 7, 7, 36, 36]
        shts = [126, 126, 126, 41, 41, 96, 96]

        # SHRG: Generalized 850-200 hPa shear magnitude (kt *10) vs time (takes
        # into account all levels from 1000 to 100 hPa
        # SHGC: Same as SHRG but with vortex removed and averaged from 0-500 km relative
        # to 850 hPa vortex center
        # should I use shrs and shts for this first test since most obs data are in the
        # lower tropsophere?

        # this value clips off problematic indices from the P-3 track for
        # the shear angle calculation. The first value is added to the lowest index, while
        # the second index is subtracted from the highest crl index
        # not really that necessary
        crl_shear_inds = ( (100, 400 ), ( 50, 50), ( 0, 0), (0, 0), (0, 0),
                ( 0, 400), ( 50, 0))

        shear_mag = shrd
        shear_dir = shtd
        shear_quads = [ ('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL'),
                ('DL', 'UR'), ('UR', 'DL'), ('DR', 'UL'), ('UL', 'DR')]

        ships_lat = [ 14.5, 14.5, 14.5, 16.5, 16.5, 20.3, 20.3]
        ships_lon = [ 50.6, 50.6, 50.6, 52.9, 52.9, 58.0, 58.0]
        ships_tlat = [ 14.4, 14.4, 14.4, 16.4, 16.4, 20.4, 20.4]
        ships_tlon = [ 50.6, 50.6, 50.6, 52.9, 52.9, 52.9, 58.0, 58.0]
        tc_center_lat = ships_tlat
        tc_center_lon = ships_tlon

        intensity = [135, 135, 135, 105, 105, 120, 120 ]
        intensity_cat = find_cat( intensity)


    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'in_situ_path': in_situ_path,
        'crl_list': crl_list, 'tdr_list': tdr_list, 'in_situ_list': in_situ_list,
        'new_crl_path': new_crl_path, 'new_tdr_path': new_tdr_path,
        'new_flight_data_path': new_flight_data_path,
        'eyewall_dists': eyewall_dists, 'in_situ_eyewall_dists': in_situ_eyewall_dists,
        'crl_range': crl_range, 'xlims': xlims, 'xtype': xtype,
        'dates': dates, 'eye_pass': eye_pass, 'tc_name': tc_name,
        'shift': shift, 'stretch': stretch,
        'shear_dir': shear_dir, 'shear_mag': shear_mag, 'tc_center_lat': tc_center_lat,
        'tc_center_lon': tc_center_lon,
        'crl_shear_inds': crl_shear_inds,
        'um_crl_path': updated_matrices_crl_path,
        'intensity': intensity, 'intensity_cat': intensity_cat,
        'shear_quads': shear_quads}

    return tcdata



# this function simply looks at the date provided in the tcdata dictionary and
# returns the appropriate name of the crl dataset
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
    elif date == '09-25':
        crl_data = crl_list[ 15]
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


def choose_in_situ_date( date, insitu_list):
    if date == '08-16':
        insitu_data = insitu_list[ 7]
    elif date == '08-17':
        insitu_data = insitu_list[ 8]
    elif date == '08-18':
        insitu_data = insitu_list[ 9]
    elif date == '08-19':
        insitu_data = insitu_list[ 12]
    elif date == '08-20':
        insitu_data = insitu_list[ 14]
    elif date == '08-21':
        insitu_data = insitu_list[ 16]
    elif date == '08-27':
        insitu_data = insitu_list[ 17]
    elif date == '08-28':
        insitu_data = insitu_list[ 18]
    elif date == '08-29':
        insitu_data = insitu_list[ 20]
    elif date == '09-26':
        insitu_data = insitu_list[ 31]
    elif date == '09-27':
        insitu_data = insitu_list[ 32]
    elif date == '09-29':
        insitu_data = insitu_list[ 34]
    else:
        print( 'update if statement!')
        return
    return insitu_data


def choose_tdr_data( tc_name, tdr_list, counter):
    if tc_name.casefold() == 'grace':
        inbound_data = [ tdr_list[ 4], tdr_list[ 6], tdr_list[ 8],
            tdr_list[ 10], tdr_list[ 12], tdr_list[ 14],
            tdr_list[ 16], tdr_list[ 18], tdr_list[ 20] ]
        outbound_data = [ tdr_list[ 5], tdr_list[ 7], tdr_list[ 9],
            tdr_list[ 11], tdr_list[ 13], tdr_list[ 15],
            tdr_list[ 17], tdr_list[ 19], tdr_list[ 21] ]
    elif tc_name.casefold() == 'henri':
        inbound_data = [ tdr_list[ 0], tdr_list[ 4], tdr_list[ 6], tdr_list[ 8], tdr_list[ 10] ]
        outbound_data = [tdr_list[ 1], tdr_list[ 5], tdr_list[ 7], tdr_list[ 9], tdr_list[ 11] ]
    elif tc_name.casefold() == 'ida':
        inbound_data = [ tdr_list[0], tdr_list[2], tdr_list[12], tdr_list[ 20]]
        outbound_data = [ tdr_list[1], tdr_list[3], tdr_list[13], tdr_list[ 21]]
    elif tc_name.casefold() == 'sam':
        inbound_data = [ tdr_list[ 0], tdr_list[ 2], tdr_list[ 4], # tdr_list[ 6],
            tdr_list[ 8], tdr_list[ 10], tdr_list[ 12], tdr_list[ 14] ]
        outbound_data = [ tdr_list[ 1], tdr_list[ 3], tdr_list[ 5], # tdr_list[ 7],
            tdr_list[9], tdr_list[ 11], tdr_list[ 13], tdr_list[ 15] ]
    else:
        print( 'update if statement!')
        return 1, 1
    return inbound_data[ counter], outbound_data[ counter]


# the function below is useful for loading the new crl and tdr datasets created
# by the scripts found in save-new-datasets!
# they can be found under C:\Users\etmu9498\research\data\tdr-new

# inputs: tc_name is the name of the storm, and counter looks at a specific eye pass
# from that storm
def choose_new_data( tc_name, counter):
    tdr_path = "/Users/etmu9498/research/data/tdr-new"
    tdr_list = make_plots.load_tdr( tdr_path, print_files=False)
    crl_path = "/Users/etmu9498/research/data/crl-new"
    crl_list = make_plots.load_crl( crl_path, print_files=False)

    if tc_name.casefold() == 'grace':
        tdr_data = [  tdr_list[ 0], tdr_list[ 0], tdr_list[ 0], # this dataset hasn't been created for 8/17 yet :(
            tdr_list[ 0], tdr_list[ 1], tdr_list[ 2],
            tdr_list[ 3], tdr_list[ 4], tdr_list[ 5] ]
        crl_data = [ crl_list[ 0], crl_list[ 1], crl_list[ 2],
            crl_list[ 3], crl_list[ 4], crl_list[ 5],
            crl_list[ 6], crl_list[ 7], crl_list[ 8]  ]
    elif tc_name.casefold() == 'henri':
        tdr_data = [ tdr_list[ 6], tdr_list[ 7], tdr_list[ 8], tdr_list[ 9], tdr_list[ 10]  ]
        crl_data = [ crl_list[ 9], crl_list[ 10], crl_list[ 11], crl_list[ 12], crl_list[ 13] ]
    elif tc_name.casefold() == 'ida':
        tdr_data = [ tdr_list[ 11], tdr_list[ 12], tdr_list[ 13], tdr_list[ 14]]
        crl_data = [ crl_list[ 14], crl_list[ 15], crl_list[ 16], crl_list[ 17]]
    elif tc_name.casefold() == 'sam':
        tdr_data = [ tdr_list[ 15], tdr_list[ 16], tdr_list[ 17], tdr_list[ 18],
            tdr_list[ 19], tdr_list[ 20], tdr_list[ 21] ]
        crl_data = [ crl_list[ 18], crl_list[ 19], crl_list[ 20], crl_list[ 21],
            crl_list[ 22], crl_list[ 23], crl_list[ 24] ]
    else:
        print( 'update if statement!')
        return 1, 1
    return tdr_data[ counter], crl_data[ counter]


def choose_new_in_situ_name( tc_name, counter):
    path = "/Users/etmu9498/research/data/in-situ-new"
    list = make_plots.load_flight_level( path, print_files=False)

    if tc_name.casefold() == 'grace':
        name = [ list[ 0], list[ 1], list[ 2],
            list[ 3], list[ 4], list[ 5],
            list[ 6], list[ 7], list[ 8] ]
    elif tc_name.casefold() == 'henri':
        name = [ list[ 9], list[ 10], list[ 11], list[ 12], list[ 13] ]
    elif tc_name.casefold() == 'ida':
        name = [ list[ 14], list[ 15], list[ 16], list[ 17]]
    elif tc_name.casefold() == 'sam':
        name = [ list[ 18], list[ 19], list[ 20], list[ 21],
            list[ 22], list[ 23], list[ 24] ]
    else:
        print( 'update if statement!')
        return 1, 1
    return name[ counter]

# split tcs into four categories: tropical depression, tropical storm, weak hurricane, and strong hurricane
# depending on their intensities from the NOAA document!
# input: a list of intensity values for this particular tc
def find_cat( intensity_list):
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
