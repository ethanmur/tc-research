
# make plot background white
plt.figure( figsize=(18, 18), facecolor='w')

# search whole directory for a phrase in atom
ctrl + shift + F

# find closest crl distances to a 'padding' variable: save their indices and corresponding distances!
idx1 = (np.abs(crl_data.in_situ_distance + padding )).argmin().values
idx2 = (np.abs(crl_data.in_situ_distance - padding )).argmin().values
rmw_lim1 = crl_data.rmw_negatives[ idx1]
rmw_lim2 = crl_data.rmw_negatives[idx2]

# find and interpolate between nans!
# taken from find_crl_distance.py
def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]
nans, x = nan_helper( speed)F
speed[ nans] = np.interp( x( nans), x( ~nans), speed[~nans])

# print all numpy array values, and reset things back to the default!
np.set_printoptions(threshold=np.inf)
np.set_printoptions(threshold=1000)

# is a value a nan or inf?
np.isnan()
np.isinf()

# split crl distances into positive and negative segments
dist_pos = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance > 0.0)].values
dist_neg = crl_data.in_situ_distance[ np.where( crl_data.in_situ_distance <= 0.0)].values

# add a box around a text() object!
bbox={'facecolor': 'w', 'alpha': 0.85, 'pad': 10},
verticalalignment='bottom', horizontalalignment='left',
transform=a0.transAxes)


# change number of x axis labels / grid lines
plt.locator_params(axis='x', nbins=11)

# make subplots different sizes!
fig, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [1, 2, .75]}, figsize=(12, 18), facecolor='w')

# hide an axis to prevent overlapping axes!
plt.grid(False)
plt.axis('off')


# add a mask to data with nans and then fit with linear regression!
mask = ~np.isnan(varx) & ~np.isnan(vary)
slope, intercept, r_value, p_value, std_err = stats.linregress(varx[mask], vary[mask])


# wrap code in these two lines to get rid of warnings!
warnings.filterwarnings("ignore")
warnings.filterwarnings("default")

# change font size of ticks for both axes!
a0.tick_params(axis='both', which='major', labelsize=smallfont)

# a nice succinct for loop!
[ str( item) for item in center_dists]

# take a slice aka truncate a colorbar! using a helper script saved locally
cmap = plt.get_cmap('Greys')
cmap = helper_fns.truncate_colormap(cmap, 0.2, 1.0)


# keep only overlapping elements from two numpy arrays
dist_inds = np.intersect1d( dist_inds1, dist_inds2)
