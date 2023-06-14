# step 1: calculate mean flight level temperature using all valid data

# find data limits: where P-3 is flying at 3 km.
fllim0 = np.nanargmin( np.abs(fl_data.time - lim0 ))
fllim1 = np.nanargmin( np.abs(fl_data.time - lim1 ))

# choose which variable to look at- Temperature here
# trim to where P-3 is only flying within the TC. Not during ferrying to and from eye
flfield = fl_data['TA.d'][fllim0:fllim1]

# get mean crl and fl data for all valid locations and add to correct list
flmean = np.nanmean(flfield)

# step 2: calculate man flight level temperature for just one eye pass

# eye limits are stored here
metadata = eyewall_metadata.all_metadata()
lims = metadata[year]['eyewall_limits'][date][eye_pass]
eyelim0, eyelim1 = lims[0], lims[1]

# find data limits: P-3 is within TC Sam's eye
fleyelim0 = np.nanargmin( np.abs(fl_data.time - eyelim0 ))
fleyelim1 = np.nanargmin( np.abs(fl_data.time - eyelim1 ))

# trim temperature data to TC Sam's eye pass
    flfield = fl_data['TA.d'][fleyelim0:fleyelim1]

# find the mean temperature within the eye    
flmean = np.nanmean(flfield)

# find and add cloud heights too
H = crl_data.height
power = crl_data.P_ch1[ crleyelim0:crleyelim1, :]
axis = crl_data.time[ crleyelim0:crleyelim1]
p3_height = crl_data.p3_height[ crleyelim0:crleyelim1]
cloudheights, cloudtime = find_cloud_tops.find_cloud_heights( H, power, axis, p3_height, cutoff_power = -30)
