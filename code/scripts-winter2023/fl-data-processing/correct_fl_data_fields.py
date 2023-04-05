# Author: Ethan Murray
# Created: 3/8/23
# Edited: 3/8/23
# Automatically process important flight level variables like temperature, wv, etc.

# import...
import os
import numpy as np
import matplotlib.pyplot as plt

# process one flight level variable based on its input name / parameters
# inputs:
# var_name: the name of the variable
# var: the actual data of the variable of interest
# p3_height: aircraft height array
# return: the processed variable data array!
def one_field( var_name, var, p3_heights, bad_inds, time, tcname, temp):

    # step 0:
    # trim every variable by the bad indices to get down to the desired height limits!
    var[ bad_inds] = np.nan

    # mixing ratio: remove anomalous wv peaks! up to 60 g/kg :/
    if var_name == 'MR.d':

        # use the error helper function to trim out the unrealistic wvmr peaks!
        error_lim = 30 # cut out wvmr values above 30 g/kg... sensible for tcs at 700 hPa
        var = error_helper( var_name, var, p3_heights, bad_inds, time, tcname, temp, error_lim)
        return var

    # repeat for equivalent potential temp, but with a higher cutoff!
    elif var_name == 'THETAE.d':
        error_lim = 375 # cut out values above 275 K.. sensible for tcs at 700 hPa
        var = error_helper( var_name, var, p3_heights, bad_inds, time, tcname, temp, error_lim)
        return var        

    # repeat for relative humidity
    elif var_name == 'HUM_REL.d':
        error_lim = 100 # cut out values above 100%
        var = error_helper( var_name, var, p3_heights, bad_inds, time, tcname, temp, error_lim)
        return var        

    # no errors for this variable! just return the height trimmed data
    else:
        return var


def error_helper( var_name, var, p3_heights, bad_inds, time, tcname, temp, error_lim):
    # step 1:
        # identify the large error regions
        # get large wv error peaks: more values are corrupted, though, so we need to find those, too!
        error_inds = np.where( var > error_lim)[0]
        # group indices into separate chunks! then only take the first and last inds
        # using 0 and -1
        error_inds_grouped = []
        for i, value in enumerate( error_inds):
            # start case
            if i == 0:
                error_inds_grouped.append( [ value])
            # continue to add to list case: ind difference is 1!
            elif value - error_inds[ i-1] == 1:
                error_inds_grouped[-1].append( value)
            # difference is greater than one case: start a new chunk!
            elif value - error_inds[ i-1] > 1:
                error_inds_grouped.append( [value])
        # go through the grouped data and save pairs of starting and ending values!
        error_pairs = []
        for i, listi in enumerate( error_inds_grouped):
            error_pairs.append( [ listi[0], listi[-1]])

        # the error pairs look surprisingly great!!
        #print( error_inds)
        #print( error_inds_grouped)
        #print( error_pairs)

        # step 2:
        # using the pair inds, search on either side of the large error regions
        # for more error values! Even if wv = 15, it still could be an unnatural, corrupted
        # peak. so, search for a trough in the data, and call that point where things go
        # back to normal :)
        # make sure to account for corner / edge cases!

        # save new pairs here
        error_pairs_updated = []
        for i in range( len( error_pairs)):
            error_pairs_updated.append( [])
        # do this for every error pair: save new, wider pairs below
        for pairi, pairval in enumerate( error_pairs):
            lefti, righti = pairval[ 0], pairval[ 1]

            # one pair code
            # search to the left of the first error data point
            newlefti = lefti
            search = True
            while search:
                # reached the first index case! just save that val
                if newlefti == 0:
                    error_pairs_updated[ pairi].append( newlefti)
                    search = False
                # break case: wv value is going up again! cut data off here
                if var[ newlefti] - var[ newlefti + 1] > 0.0:
                    error_pairs_updated[ pairi].append( newlefti)
                    search = False
                # pass case: keep searching for edge of wv errors
                else:
                    newlefti -= 1

            # search to the right of the first error data point
            newrighti = righti
            search = True
            while search:
                # reached the last index case! just save that val
                if newrighti == len( var) - 1:
                    error_pairs_updated[ pairi].append( newrighti)
                    search = False
                # break case: wv value is going up again (to the right)! cut data off here
                if var[ newrighti ] - var[ newrighti - 1] > 0.0:
                    error_pairs_updated[ pairi].append( newrighti)
                    search = False
                # pass case: keep searching for edge of wv errors
                else:
                    newrighti += 1

        # step 3:
        # make arrays stretching from one error index to another!
        error_arrays = []
        for pairi, pairval in enumerate( error_pairs_updated):
            error_arrays.append( np.arange( pairval[0], pairval[1] + 1).tolist() )
        # flatten the list of lists!
        error_list = [item for sublist in error_arrays for item in sublist]

        # step 4:
        # trim the wv data to get rid of error regions!
        # by "trim", i mean fill these data points with nans... actually shrinking the array
        # would be a pain
        # no need to trim other data: wv is often corrupted just on its own
        var[ error_list] = np.nan

        p3_heights[ error_list] = np.nan
        temp[ bad_inds] = np.nan

        return var
