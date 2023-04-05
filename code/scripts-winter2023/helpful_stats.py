import os
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import math
import seaborn as sns
import scipy

# return high and low confidence limits (and the mean) for a specified time series
# and confidence interval. Use t statistics to do so! Good for small (N<30) samples.
def t_test_intervals( data, confidence=.95, nstar=False):
    diff = 1 - confidence
    tableval = confidence + diff / 2
    data_inds = np.where( ~ np.isnan( data)) [0]
    data = np.array( data) [data_inds]

    if nstar:
        df = nstar
    else:
        N = len( data)
        df=N-1

    mean = data.mean()
    std = data.std()
    # find the t test upper and lower limits
    tstat = scipy.stats.t.ppf( tableval, df)
    low_limitt = mean - tstat * ( std / np.sqrt( N-1))
    high_limitt = mean + tstat * ( std / np.sqrt( N-1))
    return mean, low_limitt, high_limitt

# return high and low confidence limits (and the mean) for a specified time series
# and confidence interval. Use z statistics to do so! Good for large ( N>30) samples.
def z_test_intervals( data, confidence=.95):
    diff = 1 - confidence
    tableval = confidence + diff / 2
    data_inds = np.where( ~ np.isnan( data)) [0]
    data = np.array( data) [data_inds]
    N = len( data)
    df=N-1
    mean = data.mean()
    std = data.std()
    # find the z test upper and lower limits
    zstat = scipy.stats.norm.ppf( tableval)
    low_limitz = mean - zstat * ( std / np.sqrt( N))
    high_limitz = mean + zstat * ( std / np.sqrt( N))
    return mean, low_limitz, high_limitz


# autocorrelation functions- both should return the same answer
#Method #1
#Calculate the autocorrelation using numpy correlate lagN
def method1( t1_m, t2_m, N, lag, sigma):
    return np.correlate( t1_m, t2_m, mode='valid')/( N - lag)/( sigma**2)
#Method #2
#Calculate the autocorrelation using numpy dot (direct calculation)
## (https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.dot.html)
## Barnes Chapter 2 Eq. 68 divided by the variance
def method2( t1_m, t2_m, N, lag, sigma):
    return np.dot(t1_m,t2_m)/(N-lag)/sigma**2


# Nstar functions
## Calculate effective sample size (N*)
## Note: Leith function to estimate N* is not well behaved for large autocorrelations
## Prof. Kay recommends using Wilks
def nstar( data):

    # put data into correct format (xarray)
    if type( data) != type( xr.DataArray()):
        data = xr.DataArray( data)
        
    N = len( data)
    sigma = np.std( data).values  ## calculate the standard deviation
    xmean = np.mean( x_norm).values  ## calculate the mean
    lag=1
    x_norm_pandas = x_norm.to_pandas()

    x_t1_m = x_norm_pandas.iloc[0 : -1 * lag] - xmean
    x_t2_m = x_norm_pandas.iloc[ lag : ] - xmean

    # find the lag 1 autocorrelation!
    AR1 = method1( x_t1_m, x_t2_m, xN, lag, xsigma)
    # find Nstar!
    Nstar_wilks=round(((1-AR1)/(1+AR1))*N) ## Barnes Chapter 2 eq. 88
    # Nstar_leith=round((-0.5*np.log(xAR1))*xN) ## Barnes Chapter 2 eq. 90
    return Nstar_wilks
