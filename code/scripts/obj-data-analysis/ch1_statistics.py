import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import xarray as xr
import math
import seaborn as sns
import scipy

os.chdir(  "/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import make_plots_new_heights
import helper_fns
import cloud_height

os.chdir(  "/Users/etmu9498/research/code/scripts/plotting")
import simple_flight_level_plot



def t_test_intervals( data, confidence=.95):
    diff = 1 - confidence
    tableval = confidence + diff / 2

    data_inds = np.where( ~ np.isnan( data)) [0]
    data = np.array( data) [data_inds]
    N = len( data)
    df=N-1

    mean = data.mean()
    std = data.std()
    # find the t test upper and lower limits
    tstat = scipy.stats.t.ppf( tableval, df)
    low_limitt = mean - tstat * ( std / np.sqrt( N-1))
    high_limitt = mean + tstat * ( std / np.sqrt( N-1))

    return mean, low_limitt, high_limitt



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




#####################
## old methods i've been using prior to starting objective data analysis!
## they might be wrong and too forgiving lol
#####################

# confidence interval function taken from stack overflow lol
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.nanmean(a), scipy.stats.sem(a, nan_policy='omit')
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

# from the same stack overflow thread as above!
def confidence2(data, confidence=0.95):
    a = 1.0 * np.array(data)
    out = scipy.stats.t.interval( confidence, len(a)-1, loc=np.mean(a), scale=scipy.stats.sem(a))
    low_bound, high_bound = out[ 0], out[ 1]

    return np.mean( a), low_bound, high_bound

# another method for finding confidence intervals, this time taken from a different stack overflow post
# matches the two functions above
def confidence3( data, confidence=.95):
    a = 1.0 * np.array(data)
    N = len( a)
    mean, sigma = np.mean(a), np.std(a)
    out = scipy.stats.norm.interval( confidence, loc=mean, scale=sigma/math.sqrt(N))
    low_bound, high_bound = out[ 0], out[ 1]
    return mean, low_bound, high_bound

# this function finds the confidence intervals for a single draw (one random choice)
# rather than a mean draw from lots of selections. The confidence interval is way worse
# in this case :( but idk if this is the statistic that I want to use
def single_draw_confidence( data, confidence=.95):
    a = 1.0 * np.array(data)
    N = len( a)
    mean, sigma = np.mean(a), np.std(a)
    out = scipy.stats.norm.interval( confidence, loc=mean, scale=sigma)
    low_bound, high_bound = out[ 0], out[ 1]
    return mean, low_bound, high_bound
