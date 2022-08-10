# a file containing x axes, dates, and other metadata for plotting tc data automatically!
# other python scripts call these functions to load in data as dictionaries
import os
os.chdir( "/Users/etmu9498/research/code/scripts/")
import make_plots


# this function is called by make_plots.plot_full_datasets()!!
def plot_all_crl_data( tc='sam'):

    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    crl_list = make_plots.load_crl(crl_path, print_files=False)
    if tc.casefold() == 'fred':
        dates = ['08-11', '08-12-am', '08-12-pm', '08-13' ]
        tc_name = 'Fred'

    elif tc.casefold() == 'grace':
        dates = [ '08-16', '08-17', '08-18', '08-19' ]
        tc_name = 'Grace'

    elif tc.casefold() == 'henri':
        dates = [ '08-20', '08-21' ]
        tc_name = 'Henri'

    elif tc.casefold() == 'ida':
        dates = [ '08-27', '08-28', '08-29']
        tc_name = 'Ida'

    elif tc.casefold() == 'sam':
        dates = ['09-25', '09-26', '09-27', '09-29']
        tc_name = 'Sam'
    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata
    tcdata = {
        'crl_path': crl_path, 'crl_list': crl_list, 'dates': dates, 'tc_name': tc_name }
    return tcdata



# this function is called in plot_all_tc_data.plot()!!
def all_data( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here
    # paths to data
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    in_situ_path = "/Users/etmu9498/research/data/in-situ"

    new_tdr_path = "/Users/etmu9498/research/data/tdr-new"
    new_crl_path = "/Users/etmu9498/research/data/crl-new"
    # load a list of available data to these variables
    crl_list = make_plots.load_crl(crl_path, print_files=False)
    in_situ_list = make_plots.load_flight_level(in_situ_path, print_files=False)

    if tc.casefold() == 'grace':
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        crl_range = [
            (0, 2100), (2800, 3570, 3655, 4400), (5000, 6200), # 8/18
            # eye 2 had a spiral, 4 indices get rid of it
            (0, 1600), (2100, 3300), (4000, 5700) ] # 8/19 # eye 2 (2100, 3900) removed... add it again!
        xlims = [
            ( 22, 18), (-86.5, -82.5 ), (21, 18 ),
            (-92, -89), (-92, -90), (-93, -89) ] # , (-92, -90)
        xtype = [ 'lat', 'lon', 'lat', 'lon', 'lon', 'lon']

        eyewall_dists = [ ( -12.5, 17.5), (-35, 12.5  ), (-10, 35),
            (-62.5, -10), (10, 85), ( -55, 25) ]

        dates = [
            '08-18', '08-18', '08-18',
            '08-19', '08-19', '08-19' ]
        eye_pass = [ "1", "2", "3", "1", "2", "3"]

        # parameters used for getting an accurate distance axis in new crl data
        stretch = [ 1, 1, 1, 1, 1, 1]
        shift = [ 20, 0, 0, 0, 0, 10]
        tc_name = 'Grace'

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

        dates = [ '08-20', '08-20',
            '08-21', '08-21', '08-21' ]
        stretch = [ 1, 1, 1, 1, 1] # .66667
        shift = [ -50, 0, 0, 0, -5] # -12.5

        eye_pass = [ "1", "3", "1", "2", "3"]
        tc_name = 'Henri'

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

        dates = [ '08-27', '08-27', '08-27', '08-29']
        stretch = [ 1, 1, 1, 1]
        shift = [ 0, 0, 0, 90]

        eye_pass = [ "1", "2", "7", "2"]
        tc_name = 'Ida'

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

        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-29', '09-29']
        stretch = [ 1, 1, 1, 1, 1, 1, 1]
        shift = [ 65, 0, 0, 0, 0, 0, 0]

        eye_pass = [ "1", "2", "3", "2", "3", "1", "2"]
        tc_name = 'Sam'

    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'in_situ_path': in_situ_path,
        'crl_list': crl_list, 'tdr_list': tdr_list, 'in_situ_list': in_situ_list,
        'new_crl_path': new_crl_path, 'new_tdr_path': new_tdr_path, 'eyewall_dists': eyewall_dists,
        'crl_range': crl_range, 'xlims': xlims, 'xtype': xtype,
        'dates': dates, 'eye_pass': eye_pass, 'tc_name': tc_name,
        'shift': shift, 'stretch': stretch}
    return tcdata


# tc='sam' means that the function defaults to running on sam's data if no tc
# name is provided!
def choose_data_eye_cloud_heights( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here

    if tc.casefold() == 'grace':
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
            (0, 1565, 1650, 2400), (3100, 5100), (6700, 8100) ] # 8/21
            # eye 1 had a spiral, 4 indices get rid of it

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
        crl_range = [ ( 1000, 2600), (4800, 6400), (10200, 12000), (1500, 4700)]
        xlims = [ (-85, -81), (-85, -82), (-86, -82), (-92, -87)]

        dates = [ '08-27', '08-27', '08-27', '08-29']
        eye_pass = [ "1", "2", "7", "2"]
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





# tc='sam' means that the function defaults to running on sam's data if no tc
# name is provided!
def choose_data_cloud_tops_good_data( tc='sam'):
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here
    # paths to data
    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    in_situ_path = "/Users/etmu9498/research/data/in-situ"
    # load a list of available data to these variables
    crl_list = make_plots.load_crl(crl_path, print_files=False)
    in_situ_list = make_plots.load_flight_level(in_situ_path, print_files=False)

    if tc.casefold() == 'grace':
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        crl_range = [
            (0, 2100), (2800, 3570, 3655, 4400), (5000, 6200), # 8/18
            # eye 2 had a spiral, 4 indices get rid of it
            (0, 1600), (4000, 5700) ] # 8/19 # eye 2 (2100, 3900) removed
        xlims = [
            ( 22, 18), (-86.5, -82.5 ), (21, 18 ),
            (-92, -89), (-93, -89) ] # , (-92, -90)

        dates = [
            '08-18', '08-18', '08-18',
            '08-19', '08-19' ]
        eye_pass = [ "1", "2", "3", "1", "3"]
        tc_name = 'Grace'

    elif tc.casefold() == 'henri':
        tdr_path = "/Users/etmu9498/research/data/tdr/henri/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
        crl_range = [
            (0, 1565, 1650, 2400), (3100, 5100), (6700, 8100) ] # 8/21
            # eye 1 had a spiral, 4 indices get rid of it

        xlims = [
            (33.5, 39), (-73, -67.5), (-72, -68.5) ]

        dates = [
            '08-21', '08-21', '08-21' ]
        eye_pass = [ "1", "2", "3"]
        tc_name = 'Henri'

    elif tc.casefold() == 'ida':
        tdr_path = "/Users/etmu9498/research/data/tdr/ida/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        # list of touples representing the min and max indices used to clip crl data
        crl_range = [ (1500, 4700)]
        xlims = [ (-92, -87)]

        dates = [ '08-29']
        eye_pass = [ "2"]
        tc_name = 'Ida'

    elif tc.casefold() == 'sam':
        tdr_path = "/Users/etmu9498/research/data/tdr/sam/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)

        # list of touples representing the min and max indices used to clip crl data
        crl_range = [
            (0, 2000), (2350, 3800), (4400, 5600),
            (2600, 3900), (4800, 6500), # took eye 1 (0, 2000), out
            (0, 2000), (2100, 3800)]
        xlims = [
            (-52, -49), (-51.5, -49.75), (-51.75, -50),
            (18, 15), (-54, -52), # (-55, -51),
            (-59.5, -56.5), (-59.5, -56.5) ]

        dates = [
            '09-26', '09-26', '09-26', '09-27', '09-27', '09-29', '09-29']
        eye_pass = [ "1", "2", "3", "2", "3", "1", "2"]
        tc_name = 'Sam'

    else:
        tcdata = "selected TC name is not yet implemented"
        return tcdata

    tcdata = {
        'crl_path': crl_path, 'tdr_path': tdr_path, 'in_situ_path': in_situ_path,
        'crl_list': crl_list, 'tdr_list': tdr_list, 'in_situ_list': in_situ_list,
        'crl_range': crl_range, 'xlims': xlims,
        'dates': dates, 'eye_pass': eye_pass, 'tc_name': tc_name }
    return tcdata






def choose_data_eyewall_slope( tc='sam'):

    crl_path = "/Users/etmu9498/research/data/CRL_data/2021"
    crl_list = make_plots.load_crl(crl_path, print_files=False)
    # choose proper paths to data and formatting options
    # casefold() allows for any capitalization variation to work here

    if tc.casefold() == 'fred':
        tdr_path = "/Users/etmu9498/research/data/tdr/fred/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
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
        tdr_path = "/Users/etmu9498/research/data/tdr/grace/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
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
        tdr_path = "/Users/etmu9498/research/data/tdr/henri/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
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
        tdr_path = "/Users/etmu9498/research/data/tdr/ida/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
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

    elif tc.casefold() == 'sam':
        tdr_path = "/Users/etmu9498/research/data/tdr/sam/nc-files"
        tdr_list = make_plots.load_tdr(tdr_path, print_files=False)
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
        'tdr_list': tdr_list, 'xlims': xlims, 'dates': dates,
        'eye_pass': eye_pass, 'tc_name': tc_name }
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
        inbound_data = [ tdr_list[ 10], tdr_list[ 12], tdr_list[ 14],
            tdr_list[ 16], tdr_list[ 18], tdr_list[ 20] ]
        outbound_data = [ tdr_list[ 11], tdr_list[ 13], tdr_list[ 15],
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
        tdr_data = [ tdr_list[ 0], tdr_list[ 1], tdr_list[ 2],
            tdr_list[ 3], tdr_list[ 4], tdr_list[ 5] ]
        crl_data = [ crl_list[ 0], crl_list[ 1], crl_list[ 2],
            crl_list[ 3], crl_list[ 4], crl_list[ 5] ]
    elif tc_name.casefold() == 'henri':
        tdr_data = [ tdr_list[ 6], tdr_list[ 7], tdr_list[ 8], tdr_list[ 9], tdr_list[ 10]  ]
        crl_data = [ crl_list[ 6], crl_list[ 7], crl_list[ 8], crl_list[ 9], crl_list[ 10] ]
    elif tc_name.casefold() == 'ida':
        tdr_data = [ tdr_list[ 11], tdr_list[ 12], tdr_list[ 13], tdr_list[ 14]]
        crl_data = [ crl_list[ 11], crl_list[ 12], crl_list[ 13], crl_list[ 14]]
    elif tc_name.casefold() == 'sam':
        tdr_data = [ tdr_list[ 15], tdr_list[ 16], tdr_list[ 17], tdr_list[ 18],
            tdr_list[ 19], tdr_list[ 20], tdr_list[ 21] ]
        crl_data = [ crl_list[ 15], crl_list[ 16], crl_list[ 17], crl_list[ 18],
            crl_list[ 19], crl_list[ 20], crl_list[ 21] ]
    else:
        print( 'update if statement!')
        return 1, 1
    return tdr_data[ counter], crl_data[ counter]
