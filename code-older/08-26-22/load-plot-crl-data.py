# Author: Ethan Murray
# Date: 9/8/22

# functions to plot all crl cross sections.
# properties like relative humidity and temperature anomaly can be found and plotted
# in a similar way.

# path to data and dataset name
data_path = "/path/to/crl/data"
data_file = "P3_20210811H1_200000-224759.cdf"

# get data
os.chdir( data_path)
crl_data = xr.open_dataset( data_file)

# indices to cut off data: no data is trimmed in this example, though.
# trimming data can save time while plotting and prevent overlapping data.
i1 = 0
i2 = len( crl_data.time)

# which variable should be used on the x axis?
xaxis_name = 'time'

# helper function to choose x axis data
def x_axis_helper( crl_data, iex1, i2, xaxis):

    if xaxis == 'lon':
        xaxis = crl_data.Lon[i1:i2]
        return xaxis

    elif xaxis == 'lat':
        xaxis = crl_data.Lat[i1:i2]
        return xaxis

    elif xaxis == 'time':
        xaxis = crl_data.time[i1:i2]
        return xaxis


# plot CRL temperature data
def plot_T( crl_data, i1, i2, xaxis_name):
    # block out some annoying error messages
    warnings.filterwarnings("ignore")

    # choose a colormap for prettier plots
    color_map = plt.cm.get_cmap( "RdYlBu").reversed()

    # looking at different colormaps
    # color_map = plt.cm.get_cmap( "RdYlBu").reversed()
    # color_map = 'viridis'
    # color_map = plt.cm.get_cmap( "Spectral").reversed()

    # calculate temperature
    # only look at temperatures below 50 degrees C; anything above this is an error
    # .transpose() flips T matrix so that plotting with pcolormesh() works
    temp = crl_data.T[i1:i2, :].where( crl_data.T[i1:i2, :].values < 50).transpose()

    # choose an x axis
    xaxis = x_axis_helper( crl_data, i1, i2, xaxis_name)

    # plot and make things pretty
    # colorbar is from 5 C to 35 C here
    # also setting a special colormap here!
    plt.pcolormesh( xaxis, - crl_data.H, temp, cmap = color_map, vmin=5, vmax=35 )
    plt.ylabel( 'Height (km)')
    plt.xlabel( 'Time (Hours, UTC)')

    # add black background when no data is presnt
    ax = plt.gca()
    ax.set_facecolor('k')

    plt.colorbar(label="Temperature ( C)")
    plt.grid( True)

    # return to default warnings
    warnings.filterwarnings("default")

# plot water vapor data
def plot_wvmr( crl_data, i1, i2, xaxis_name):
    # choose x axis with helper script
    xaxis = x_axis_helper( crl_data, i1, i2, xaxis_name)

    # calculate wvmr
    # remove 0 and values over 20 g/kg; maybe increase this upper limit
    step1 = crl_data.WVMR.where( crl_data.WVMR.values != 0)
    step2 = step1.where( step1.values < 20)
    # trim and transpose data
    crl_wvmr = step2[i1:i2, :].transpose()

    # plot
    plt.pcolormesh( xaxis, - crl_data.H, crl_wvmr, vmin = 0, vmax =20)
    plt.colorbar(label="WVMR ( g/kg)")
    plt.ylabel( 'Height (km)')
    plt.xlabel( 'Time (Hours, UTC)')

    # plt.xlim( [ 20, 24] ) # example x data limit
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')

# plot light scattering ratio
def plot_lsr( crl_data, i1, i2, xaxis_name):
    xaxis = x_axis_helper( crl_data, i1, i2, xaxis_name)

    # calculate lsr
    step1 = crl_data.LSR[i1:i2, :].where( crl_data.LSR[i1:i2].values < 10).transpose()
    crl_lsr = step1.where( step1.values > .1)

    # plot
    plt.pcolormesh( xaxis, - crl_data.H, crl_lsr)
    plt.colorbar(label="LSR ( unitless)")
    plt.ylabel( 'Height (km)')

    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')

# plot backscattered power to channel 1
def plot_power_ch1( crl_data, i1, i2, xaxis_name, cutoff=-30):
    xaxis = x_axis_helper( crl_data, i1, i2, xaxis_name)

    # calculate power backscattered to channel 1
    # take the log of values (put into dbz)
    step1 = 10 * np.log10( crl_data.P_ch1 )
    # remove values below -30 dBz: not necessary, but backscatter under -30 dBz
    # is annoying for cloud height calculations, etc, so I usually remove it
    step2 = step1.where( step1.values > cutoff)
    crl_pch1 = step2[i1:i2, :].transpose()

    plt.pcolormesh(  xaxis, - crl_data.H, crl_pch1, vmin = cutoff, vmax =-10)
    if show_colorbar:
        plt.colorbar(label="Backscattered Ch 1 power ( dBz)")
    plt.ylabel( 'Height (km)')
    plt.grid( 'on')
    ax = plt.gca()
    ax.set_facecolor('k')
