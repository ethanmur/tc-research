# tc_metadata.py
# edited 12/16/22
# a file containing x axes, dates, and other metadata for working with tc data automatically!
# the top scripts choose and load in specific datasets as dictionaries
import os
os.chdir( "/Users/etmu9498/research/code/scripts/")
import make_plots


# this function simply looks at the date provided in the metadata dictionary and
# returns the appropriate name of a given crl dataset.
# this provides the original crl dataset names; see the choose_new_data()
# function to load new crl names
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


# this function works like the one above, but loads filenames for in situ data
def choose_in_situ_date( date, insitu_list):
    if date == '08-11':
        insitu_data = insitu_list[ 3]
    elif date == '08-12-am':
        insitu_data = insitu_list[ 4]
    elif date == '08-12-pm':
        insitu_data = insitu_list[ 5]
    elif date == '08-13':
        insitu_data = insitu_list[ 6]
    elif date == '08-16':
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


# this function looks at a tc name, a list of all tdr datasets for that name,
# and a counter to pick and return an inbound and outbound dataset.
# The code for TC fred is new, easier to update, and more efficient; this code
# should be applied to all other cases if possible!
def choose_tdr_data( tc_name, tdr_list, counter):

    if tc_name.casefold() == 'fred':
        inbound_data = []
        outbound_data = []

        # manually choose i range to pick the correct tdr datasets
        # 18 datasets, so 18 was chosen to be inclusive
        # step = 2 spaces out the inbound / outbound datasets correctly!
        # start on 4 to skip 8/11 data, not used in this analysis
        for i in range( 4, 18, 2):
            inbound_data.append( tdr_list[ i])
            outbound_data.append( tdr_list[ i + 1])

    # old way of doing things (inefficient)
    # since each tc has its own folder, the code below was kept the same :)
    elif tc_name.casefold() == 'grace':
        inbound_data = [ tdr_list[ 0], tdr_list[ 2],
            tdr_list[ 4], tdr_list[ 6], tdr_list[ 8],
            tdr_list[ 10], tdr_list[ 12], tdr_list[ 14],
            tdr_list[ 16], tdr_list[ 18], tdr_list[ 20] ]
        outbound_data = [ tdr_list[ 1], tdr_list[ 3],
            tdr_list[ 5], tdr_list[ 7], tdr_list[ 9],
            tdr_list[ 11], tdr_list[ 13], tdr_list[ 15],
            tdr_list[ 17], tdr_list[ 19], tdr_list[ 21] ]
    elif tc_name.casefold() == 'henri':
        inbound_data = [ tdr_list[ 0], tdr_list[ 4], tdr_list[ 6], tdr_list[ 8], tdr_list[ 10] ]
        outbound_data = [tdr_list[ 1], tdr_list[ 5], tdr_list[ 7], tdr_list[ 9], tdr_list[ 11] ]
    elif tc_name.casefold() == 'ida':
        inbound_data = [ tdr_list[0], tdr_list[2], tdr_list[4], tdr_list[12], tdr_list[ 20]]
        outbound_data = [ tdr_list[1], tdr_list[3], tdr_list[5], tdr_list[13], tdr_list[ 21]]
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
# they can be found under C:\Users\etmu9498\research\data\tdr-new or crl-new
# inputs: tc_name is the name of the storm, and counter looks at a specific eye pass
# from that storm.
# this code was recently edited (12/15/22) to loop through datasets instead of
# manually choosing each number. The old code snippets can be found in tc_metadata_old.py
def choose_new_data( tc_name, counter):
    tdr_path = "/Users/etmu9498/research/data/tdr-new"
    tdr_list = make_plots.load_tdr( tdr_path, print_files=False)
    crl_path = "/Users/etmu9498/research/data/crl-new"
    crl_list = make_plots.load_crl( crl_path, print_files=False)

    tdr_data = []
    crl_data = []
    if tc_name.casefold() == 'fred':
        # the for loop cycles through elements in the list of new crl
        # dataset names.
        # the range of the cycle needs to be manually set, depending on the
        # number of datasets present in the new data folders
        for i in range( 7):
            crl_data.append( crl_list[ i])
            tdr_data.append( tdr_list[ i])
    elif tc_name.casefold() == 'grace':
        for i in range( 7, 18):
            tdr_data.append( tdr_list[ i])
            crl_data.append( crl_list[ i])
    elif tc_name.casefold() == 'henri':
        for i in range( 18, 23):
            tdr_data.append( tdr_list[ i])
            crl_data.append( crl_list[ i])
    elif tc_name.casefold() == 'ida':
        for i in range( 23, 28):
            tdr_data.append( tdr_list[ i])
            crl_data.append( crl_list[ i])
    elif tc_name.casefold() == 'sam':
        for i in range( 28, 35):
            tdr_data.append( tdr_list[ i])
            crl_data.append( crl_list[ i])
    else:
        print( 'update if statement!')
        return 1, 1
    return tdr_data[ counter], crl_data[ counter]


# This code works like the function above, except that it looks at in situ data.
# old code snippets are also located in tc_metadata_old.py
def choose_new_in_situ_name( tc_name, counter):
    path = "/Users/etmu9498/research/data/in-situ-new"
    list = make_plots.load_flight_level( path, print_files=False)
    data = []
    if tc_name.casefold() == 'fred':
        for i in range( 7):
            data.append( list[ i])
    elif tc_name.casefold() == 'grace':
        for i in range( 7, 18):
            data.append( list[ i])
    elif tc_name.casefold() == 'henri':
        for i in range( 18, 23):
            data.append( list[ i])
    elif tc_name.casefold() == 'ida':
        for i in range( 23, 28):
            data.append( list[ i])
    elif tc_name.casefold() == 'sam':
        for i in range( 28, 35):
            data.append( list[ i])
    else:
        print( 'update if statement!')
        return 1, 1
    return data[ counter]


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

        crl_range = [ (0, 1650), (1650, 3300), (4200, 5500),
                (0, 1900), (2900, 4500),
                (0, 1750), (2300, 4000)]

        xlims = [
            ( -77, -72), ( -77, -72), (24, 18),
            ( -77, -73.5), (-77, -73.5 ),
            (-79, -75), (-79, -75) ]
        xtype = ['lon', 'lon', 'lat',
            'lon', 'lon',
            'lon', 'lon' ]

        # placeholders that aren't totally relevant, need to initialize variables to
        # avoid errors
        eyewall_dists = []
        stretch = []
        shift = []
        crl_shear_inds = []
        tc_center_lat = []
        tc_center_lon = []

        # older, more inclusive eyewall definitions
        in_situ_eyewall_dists = [ ( -65, 60), ( -70, 0), ( -55, 50),
            ( -45, 47.5), ( -57.5, -5),
            ( -60, 50), ( -55, 70) ]

        eyewall_dists_no_eyewalls = [ ( -10, 40), ( -70, -16), ( -55, 50),
            ( -45, 47.5), ( -52.5, -11),
            ( -57.5, 40), ( -55, 50) ]

        dates = [ '08-12-am', '08-12-am', '08-12-am',
            '08-12-pm', '08-12-pm',
            '08-13', '08-13' ]

        eye_pass = [ "1", "2", "3", "1", "2", "1", "2"]

        tc_name = 'Fred'

        # Ships data used:
        # 12 UTC on 08-12 for 08-12-am
        # 0 UTC on 08-13 for 08-12-pm
        # 12 UTC on 08-13 for 08-13
        shrd = [ 187, 187, 187, 181, 181, 179, 179]
        shtd = [ 109, 109, 109, 91, 91, 97, 97]
        shear_dir = shtd
        shear_mag = shrd

        # found on 1/4/23. determined with shear_angles.py scripts! saved here for convenience
        shear_quads = [('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL'),
                    ('UL', 'DR'), ('DL', 'UR'), ('UL', 'DR'), ('UL', 'DR')]


        intensity = [ 30, 30, 30, 32.5, 32.5, 30, 30 ]
        intensity_cat = find_cat( intensity)


    elif tc.casefold() == 'grace':
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        crl_range = [ ( 3000, 4700), ( 5000, 6658), # 8/16
            ( 0, 1800), (2250, 4100), (4800, 6550), # 8/17
            (0, 2100), (2800, 3570, 3655, 4400), (5000, 6200), # 8/18
            # eye 2 had a spiral, 4 indices get rid of it
            (0, 1600), (2100, 3300), (4000, 5700) ] # 8/19 # eye 2 (2100, 3900) removed... add it again!

        xlims = [ (-73, -69), (-72, -69),
            (20, 16), (-78.5, -74), (20, 17), # 8/17
            ( 22, 18), (-86.5, -82.5 ), (21, 18 ), # 8/18
            (-92, -89), (-92, -90), (-93, -89) ] # 8/19
        xtype = ['lon', 'lon', 'lat', 'lon', 'lat', 'lat', 'lon', 'lat', 'lon', 'lon', 'lon']

        # original eyewall distances when using tdr axis
        eyewall_dists = [ ( -12.5, 17.5), (-35, 12.5  ), (-10, 35),
            (-62.5, -10), (10, 85), ( -55, 25) ]


        # new eyewall distances for height corrected data with an x axis based on in situ data
        # new, precise limits
        in_situ_eyewall_dists = [ ( -55, -27.5), (25, 55), # temp values for 8/16
                (-65, 72.5), (11, 57), (-12.5, 53), # old pass 2 and 3 eyewalls: (-75, 57.5), (-48, 53),
                (-11.5, 17), (-22.0, 7.5), (-12.5, 26),
                (-32.5, 11), ( 10, 70), (-55.5, 29.5)]

        # the first bit of the sloping eyewall is removed, but once there's a sharp drop,
        # that's considered a part of the eye
        # the same limits as above, but with sloping eyewalls removed from the datasets!
        eyewall_dists_no_eyewalls = [ (-53, -36), (25, 55),
                            (-65, 70), (14, 54), (-11, 48), # passes 2 and 3 are tough... look at again
                            (-9, 16), (-12.5, 6.5), (-9, 24),
                            (-31, 11), (10, 68), (-52, 27)]

        # old no eyewall limits, cutting off whole eyewall... too aggressive
        # eyewall_dists_no_eyewalls = [ (-53, -36), (25, 55),
        #                     (-65, 63), (14, 52), (-7, 48),
        #                     (-9, 16), (-12, 6), (-7, 16),
        #                     (-31, 10), (12.5, 68), (-40, 27)]

        dates = [
            '08-16', '08-16',
            '08-17', '08-17', '08-17',
            '08-18', '08-18', '08-18',
            '08-19', '08-19', '08-19' ]
        eye_pass = [ "1", "2", "1", "2", "3", "1", "2", "3", "1", "2", "3"]

        # parameters used for getting an accurate distance axis in new crl data
        stretch = [ 1, 1, 1, 1, 1, 1]
        shift = [ 20, 0, 0, 0, 0, 10]

        tc_name = 'Grace'

        # Ships data used:
        # 12 UTC on 08-16 for 08-16
        # 12 UTC on 08-17 for 08-17
        # 0 UTC on 08-19 for 08-18
        # 0 UTC on 08-20 for 08-19
        shrd = [ 160, 160, 107, 107, 107, 105, 105, 105, 81, 81, 81]
        shtd = [ 154, 154, 142, 142, 142, 135, 135, 135, 169, 169, 169]

        # newly updated shear quads: 1/4/23
        shear_quads = [('UR', 'DL'), ('UL', 'DR'),
                    ('UL', 'DR'), ('DL', 'UR'), ('DL', 'UR'),
                    ('UR', 'DL'), ('DL', 'UR'), ('DR', 'UL'),
                    ('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL')]

        # shear_quads = [ (0, 0), (0, 0), # need to update these values
        #             ('UL', "DR"), ('DL', 'UR'), ('DL', 'UR'),
        #             ('UR', 'DL'), ('DL', 'UR'), ('DR', 'UL'),
        #             ('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL')] # determined with shear_angles.py scripts! saved here for convenience

        shear_dir = shtd
        shear_mag = shrd

        crl_shear_inds = []
        tc_center_lat = []
        tc_center_lon = []

        # tc intensity (kts) at time of flight, taken from noaa summary docs
        intensity = [ 32.5, 32.5, 45, 45, 45, 70, 70, 70, 55, 55, 55 ]
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

        # the same limits as above, but with sloping eyewalls removed from the datasets!
        eyewall_dists_no_eyewalls = [ (5, 45), ( 5, 36),
                            (-3, 12.5), (-12.5, 21), (-15, 32)]

        # old no eyewall limits... keep these values!
        # eyewall_dists_no_eyewalls = [ (5, 45), ( 5, 36),
        #                     (-3, 12.5), (-7.5, 10.5), (-15, 24)]


        dates = [ '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        stretch = [ 1, 1, 1, 1, 1] # .66667
        shift = [ -50, 0, 0, 0, -5] # -12.5

        eye_pass = [ "1", "3", "1", "2", "3"]
        tc_name = 'Henri'

        # Ships data used:
        # 0 UTC on 08-21 for 08-20
        # 0 UTC on 08-22 for 08-21
        shrd = [ 128, 128, 60, 60, 60]
        shtd = [ 179, 179, 240, 240, 240]
        shear_dir = shtd
        shear_mag = shrd

        # checked on 1/4/23: looks good!
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
        crl_range = [ ( 1000, 2600), (4800, 6400), (6900, 8000), (10200, 12000), # fourth eye from 8/27 has been added!
                    (1500, 4700)]
        xlims = [ (-85, -81), (-85, -82), (21, 23.5), (-86, -82), (-92, -87)]
        xtype = [ 'lon', 'lon', 'lat', 'lon', 'lon']

        # maybe add data from early 8/28? Looks like the P-3 was sampling something
        # else, but worth a closer look... maybe 1 good eye dataset here!

        eyewall_dists = [(-15, 65), (-35, 15), (-100, 100), (-30, 50), (-12.5, 10)] # fix eyewall dists!
        in_situ_eyewall_dists = [(-25, 60), (-26, 34), (-20, 25), (-15, 75), (-7, 20)]

        # the same limits as above, but with sloping eyewalls removed from the datasets!
        eyewall_dists_no_eyewalls = [ (-12.5, 57.5), (-20, 31), (-7, 19), (7, 63),
                             (-2.5, 15)]

        # old no eyewall limits... keep these values!
        # eyewall_dists_no_eyewalls = [ (-11, 57.5), (-14, 20), (-7, 8), (8, 63),
        #                      (-1, 15)]

        dates = [ '08-27', '08-27', '08-27', '08-27', '08-29']
        stretch = [ 1, 1, 1, 1]
        shift = [ 0, 0, 0, 90]

        eye_pass = [ "1", "2", "3", "7", "1"]
        tc_name = 'Ida'

        # Ships data used:
        # 0 UTC on 08-28 for 08-27
        # 18 UTC on 08-29 for 08-29
        shrd = [ 115, 115, 115, 115, 112]
        shtd = [ 61, 61, 61, 61, 127]

        crl_shear_inds = []

        shear_dir = shtd
        shear_mag = shrd

        # shear quads updated on 1/4/23
        shear_quads = [('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL'), ('DL', 'UR'), ('DL', 'UR')]


        # shear_quads = [ ('UL', 'DR'), ('DL', 'UR'), (0, 0), ('DL', 'UR'), ('DL', 'UR')]

        tc_center_lat = []
        tc_center_lon = []

        intensity = [70, 70, 70, 70, 130 ]
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

        # the same limits as above, but with sloping eyewalls removed from the datasets!
        eyewall_dists_no_eyewalls = [ (0, 29), (-1, 22), (-36, -2.5),
                            (0, 31), (-15, 5),
                            (-23, 24), (-17, 25)]

        # old no eyewall limits... keep these values!
        # eyewall_dists_no_eyewalls = [ (0, 26), (-1, 22), (-36, -6),
        #                     (0, 27.5), (-13, 2.5),
        #                     (-23, 24), (-16, 24)]

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

        # Ships data used:
        # 0 UTC on 09-28 for 09-27
        # 0 UTC on 09-29 for 09-28
        # 0 UTC on 09-30 for 09-29
        shrd = [ 85, 85, 85, 110, 110, 103, 103]
        shtd = [ 44, 44, 44, 74, 74, 57, 57]

        # old but very wrong values!!! 0 UTC for the previous day was used for all cases, not the current day
        # shrd = [ 90, 90, 90, 85, 85, 103, 103 ]
        # shtd = [ 64, 64, 64, 44, 44, 57, 57]

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


        # fixed shear quads found with correct shrd and shtd data
        # only one quadrant determination for 9/28 pass 3 changed:
        # it was within 16 degrees of the crossover point!
        shear_quads = [ [('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL'),
            ('DL', 'UR'), ('DR', 'UL'), ('DR', 'UL'), ('UL', 'DR')]]

        # old shear quads determined by too early shrd data, described above
        # shear_quads = [ ('UL', 'DR'), ('DL', 'UR'), ('DR', 'UL'),
        #         ('DL', 'UR'), ('UR', 'DL'), ('DR', 'UL'), ('UL', 'DR')]

        ships_lat = [ 14.5, 14.5, 14.5, 16.5, 16.5, 20.3, 20.3]
        ships_lon = [ 50.6, 50.6, 50.6, 52.9, 52.9, 58.0, 58.0]
        ships_tlat = [ 14.4, 14.4, 14.4, 16.4, 16.4, 20.4, 20.4]
        ships_tlon = [ 50.6, 50.6, 50.6, 52.9, 52.9, 52.9, 58.0, 58.0]
        tc_center_lat = ships_tlat
        tc_center_lon = ships_tlon

        intensity = [135, 135, 135, 105, 105, 115, 115 ]
        intensity_cat = find_cat( intensity)

    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    # add metadata to a dictionary for a given tc case!
    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'in_situ_path': in_situ_path,
        'crl_list': crl_list, 'tdr_list': tdr_list, 'in_situ_list': in_situ_list,
        'new_crl_path': new_crl_path, 'new_tdr_path': new_tdr_path,
        'new_flight_data_path': new_flight_data_path,
        'eyewall_dists': eyewall_dists, 'in_situ_eyewall_dists': in_situ_eyewall_dists, 'eyewall_dists_no_eyewalls': eyewall_dists_no_eyewalls,
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




# this separate dictionary creation script is used in plot_all_in_situ_one_day.py
# it's able to print earlier grace data and all the fred data!
# it needs a lot fewer of the fields specified in the all_data script above
def all_in_situ_metadata( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here
    # paths to data
    # load a list of available data to these variables
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    crl_list = make_plots.load_crl(crl_path, print_files=False)

    if tc.casefold() == 'fred':
        dates = [
            '08-12-am', '08-12-am', '08-12-am',
            '08-12-pm', '08-12-pm',
            '08-13', '08-13' ]
        tc_name = 'Fred'
    elif tc.casefold() == 'grace':
        dates = [
            '08-16', '08-16', '08-16',
            '08-17', '08-17', '08-17',
            '08-18', '08-18', '08-18',
            '08-19', '08-19', '08-19' ]
        tc_name = 'Grace'
    elif tc.casefold() == 'henri':
        dates = [ '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        tc_name = 'Henri'
    elif tc.casefold() == 'ida':
        dates = [ '08-27', '08-27', '08-27', '08-29']
        tc_name = 'Ida'
    elif tc.casefold() == 'sam':
        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-29', '09-29']
        tc_name = 'Sam'
    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata
    tcdata = {
        'crl_list': crl_list, 'dates': dates, 'tc_name': tc_name}
    return tcdata
