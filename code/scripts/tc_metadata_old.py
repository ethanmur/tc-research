# created 8/31/22
# a bunch of code from tc_metadata.py that hasn't been used in a while.
# I decided to throw all the extra bits into here to make the important file cleaner!
# Just put things back if tc_metadata.py isn't working lol


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
