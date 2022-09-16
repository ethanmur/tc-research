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
nans, x = nan_helper( speed)
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
