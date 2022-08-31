import os
import xarray as xr
import numpy as np

def in_situ_helper( crl_path, crl_name, cutoff_indices, return_var, float_time):

    # old way to convert strings to floats
    # faster but there's an error with empty strings
    # return_var = [ float( line) for line in return_var]

    # get rid of empty strings in height dataset
    # slower than code above but doesn't crash when there's an empty string
    return_var_temp = np.zeros( len( return_var))
    for line_ind in range( len( return_var)):
        if return_var[ line_ind] == '':
            return_var_temp[line_ind] = np.nan
        else:
            return_var_temp[ line_ind] = float( return_var[ line_ind])

    # print( len( return_var))
    # print( len( return_var_temp))
    # print( type( return_var_temp))

    return_var = return_var_temp.tolist()

    # load crl data to find the times corresponding to i1 and i2
    os.chdir( crl_path)
    crl_data = xr.open_dataset( crl_name)

    # deal with annoying indices :(
    if len( cutoff_indices) == 2:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
        time1 = crl_data.time[ i1]
        time2 = crl_data.time[ i2]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()

        return return_var[ idx1.values : idx2.values]

    elif len( cutoff_indices) == 4:
        i1 = cutoff_indices[0]
        i2 = cutoff_indices[1]
        i3 = cutoff_indices[2]
        i4 = cutoff_indices[3]
        time1 = crl_data.time[ i1]
        time2 = crl_data.time[ i2]
        time3 = crl_data.time[ i3]
        time4 = crl_data.time[ i4]

        # find the in situ times nearest the crl times
        idx1 = (np.abs(float_time - time1)).argmin()
        idx2 = (np.abs(float_time - time2)).argmin()
        idx3 = (np.abs(float_time - time3)).argmin()
        idx4 = (np.abs(float_time - time4)).argmin()

        # use the indices for nearest time values to trim down lat and lon values
        # this prevents data overlap and makes plotting faster!
        return return_var[ idx1.values : idx2.values] + return_var[ idx3.values : idx4.values]
    else:
        print( 'Error in number of indices! Update them in tc_metadata.py')
