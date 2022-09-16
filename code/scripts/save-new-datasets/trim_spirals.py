import os
import xarray as xr
import numpy as np
os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata



# this function looks through all crl datasets; if the dataset contains a microphysical
# spiral (four crl_range indices, not two), it trims it out and saves the dataset in a
# special folder
def trim_spirals( ):
    tcname_list = ['grace', 'henri', 'ida', 'sam']
    for tcname in tcname_list:
        # load the metadata for this tc
        metadata = tc_metadata.all_data( tc= tcname)
        # look at one specific cross section from this tc
        for dataset in range( len( metadata['dates'])):

            # is there a microphysical spiral?
            if len( metadata[ 'crl_range'][dataset]) == 4:
                # load the data names for this case
                crl_name = tc_metadata.choose_crl_date( metadata[ 'dates'][ dataset], metadata[ 'crl_list'])
                # load the actual crl data for editing
                os.chdir( metadata[ 'crl_path'])
                crl_data = xr.open_dataset( crl_name)

                # define i1 through i4 locally for convenience
                i1 = metadata[ 'crl_range'][dataset][0]
                i2 = metadata[ 'crl_range'][dataset][1]
                i3 = metadata[ 'crl_range'][dataset][2]
                i4 = metadata[ 'crl_range'][dataset][3]


                # get rid of the spiral, and trim each dataset down to its pass limits!
                # most of this code has been copied over from the save_crl_data.py file

                # make a copy of the crl data: this will eventually become the new, full dataset!
                crl_new = crl_data.copy()
                # drop all the coords, dims, and vars not listed!!
                crl_new = crl_new[ ['ProductionDateTime', 'VersionID']]

                # add a new attribute explaining how this dataset has been edited
                disclaimer = 'A microphysical spiral present in the middle of this dataset has been automatically removed.'
                crl_new.attrs[ 'global_att6'] = disclaimer

                # add time and height as coordinates, not variables!
                cliptime = np.concatenate( (crl_data.time[ i1 : i2 ].to_numpy(), crl_data.time[ i3:i4].to_numpy()) )
                crl_new = crl_new.assign_coords( {'time': cliptime })
                crl_new = crl_new.assign_coords( {'H': crl_data.H })

                # add time and height as coordinates, not variables!
                disttime = crl_data.time[ i1 : i4 - (i3-i2) ]
                crl_new = crl_new.assign_coords( {'time_distance_axis': disttime })
                crl_new[ 'time_distance_axis'].attrs[ 'Note:'] = "Only use this time axis for distance from TC center calculations."


                # loop through all the variables in crl_data, cropping them to the smaller i1 and i2 values provided
                var_list = list( crl_data.keys() )
                for key_ind in range( len( var_list)):

                    # saving 2D arrays
                    try:
                        # clip the data
                        currentdata = crl_data[ var_list[ key_ind]]

                        # put everything into numpy form using .data to avoid xr error
                        # .data and .to_numpy() do the same thing!
                        dataclip = np.concatenate( (currentdata[ i1 : i2 ].to_numpy(), currentdata[ i3:i4].data) )
                        # add data to the new, smaller xarray file
                        crl_new = crl_new.assign( { var_list[ key_ind]: ( ('time', 'H'), dataclip) })

                        # add all the attribute metadata to each variable in the new dataset!
                        # get the old metadata
                        attr_dict = crl_data[ var_list[ key_ind]].attrs
                        # recursively add all keys to new variable!
                        for key2 in attr_dict:
                            crl_new[ var_list[ key_ind]].attrs[ key2] = attr_dict[ key2]

                    # not a 2D array: do nothing
                    # all 1D arrays and values have been accounted for here or below
                    except (IndexError, ValueError):
                        continue

                # add lat and lon manually because it's being a pain to do it automatically
                # add them as coordinates, not variables
                lonclip = np.concatenate( (crl_data.Lon[ i1:i2].to_numpy(), crl_data.Lon[ i3:i4].to_numpy()) )
                crl_new = crl_new.assign_coords( { 'Lon':  lonclip })
                latclip = np.concatenate( (crl_data.Lat[ i1:i2].to_numpy(), crl_data.Lat[ i3:i4].to_numpy()) )
                crl_new = crl_new.assign_coords( { 'Lat':  latclip })
                # add metadata for lat and lon
                for key in ['Lon', 'Lat']:
                    attr_dict = crl_data[ key].attrs
                    # recursively add all keys to new dataset!
                    for key2 in attr_dict:
                        crl_new[ key].attrs[ key2] = attr_dict[ key2]


                # make a name for the dataset
                # filename = 'crl-' + metadata[ 'tc_name'].lower() + '-' + metadata[ 'dates'][ dataset] + '-eye-' + metadata[ 'eye_pass'][ dataset] + '.nc'
                filename = crl_name
                # save the newly created crl dataset
                crl_new.to_netcdf('/Users/etmu9498/research/data/crl-spiral-cases/' + filename)

                print( "Microphysical Spiral Removed from: " + filename)
    return
