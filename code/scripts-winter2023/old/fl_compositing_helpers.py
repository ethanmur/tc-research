# import...
import os
os.chdir("/Users/etmu9498/research/code/scripts")
import make_plots


# given a certain tc specification, this code returns lists representing the valid tc years and names
# imputs:
###################
# tc-
def get_filenames( tc):
    fl_path_root = "/Users/etmu9498/research/data/in-situ-noaa-full/"
    data_path = "/Users/etmu9498/research/data/"

    year_list = []
    tcname_list = []

    # which datasets should we make time series plots for?
    # all datasets case
    if tc == 'all':

        # make a year list (pretty standard / easy)
        os.chdir( data_path)
        year_list = [name for name in os.listdir('in-situ-noaa-full')
            if os.path.isdir(os.path.join('in-situ-noaa-full', name))]

        print( year_list)

        fl_list_count = 0

        # go through every year-
        # find the folders within that year
        # then, count the number of datasets per folder!
        os.chdir( fl_path_root)
        for yeari, yearval in enumerate( year_list):
            tcname_listi = [name for name in os.listdir( yearval)
                if os.path.isdir(os.path.join( yearval, name))]

            tcname_list.append( [])

            for namei, nameval in enumerate( tcname_listi):

                # add this name to the name list!
                tcname_list[ yeari].append( nameval)

                # go to new directory and get the count + names of files there!
                fl_listi = make_plots.load_flight_level( fl_path_root + year_list[ yeari] + '/' + nameval, print_files=False)
                fl_list_count += len( fl_listi)
        #############

        print( 'Total Number of plots to be created: ' + str(fl_list_count)+ '\n')

        return year_list, tcname_list


    # one year case- check the first two characters: a proper year input?
    elif len( tc) == 4 and ( tc[0:2] == '19' or tc[0:2] == '20'):

        # check if this year has a folder
        folder_present = False
        for folder in year_list:
            # the folder exists
            if folder == tc:
                folder_input = [ tc]
                folder_present = True

                fl_list_count = 0
                # count up all the files in the name folders
                tcname_list = [name for name in os.listdir( yearval)
                    if os.path.isdir(os.path.join( yearval, name))]
                for namei, nameval in enumerate( tcname_list):
                    # go to new directory and get the count + names of files there!
                    fl_listi = make_plots.load_flight_level( fl_path_root + tc + '/' + tcname_list[ namei], print_files=False)
                    fl_list_count += len( fl_listi)

                print( 'Total Number of plots to be created: ' + str( fl_list_count) + '\n')

        # no folder present case:
        if folder_present == False:
            print( "No folder present. Please download flight level data and add a new folder for it!")
            return

    # input names of specific runs as a dict! each year has its own entry
    elif type( tc) == type( {}):
        # make the year_list: just the dates saved in the dictionary!
        folder_input = list( tc.keys())

    # made it out of the loop?
    else:
        print( 'implement the else case! cannot handle individual plots yet.')


    return year_list, tcname_list
