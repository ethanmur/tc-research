
# load csv formatted data
# use os to change to filepath, then load by just typing filename
pd.read_csv('filename.csv')

# load and save pandas dataframes using pickles
pd.read_pickle("/Users/etmu9498/research/data/simple_distances.pkl") # load
# save
os.chdir("to/file/")
fl_df_total.to_pickle("all_passes_" + str( fl_df_total.shape[0]) + "_cases.pkl")

# change print settings for arrays
np.set_printoptions(threshold=np.inf)
np.set_printoptions(threshold=1000)

# print out a full pandas dataframe
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df)


# make unique colors for each pass / dataset
cmap = matplotlib.cm.get_cmap('RdYlBu')
evenly_spaced = np.linspace(0, 1, rmw_case_count).tolist()
rgblist = cmap(evenly_spaced)

# fancy legend options
leg = plt.legend( loc='upper right', bbox_to_anchor=(1.011, 1.32), fancybox=False, shadow=False, fontsize=18, facecolor='w', framealpha=1)

# colorbar with zero set not at midpoint
from matplotlib import colors
divnorm=colors.TwoSlopeNorm(vmin=-5., vcenter=0., vmax=10)
pcolormesh(your_data, cmap="coolwarm", norm=divnorm)