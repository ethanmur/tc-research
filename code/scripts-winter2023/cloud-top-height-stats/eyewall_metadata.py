import os
import numpy as np
import xarray as xr

os.chdir("/Users/etmu9498/research/code/scripts")
import tc_metadata
import make_plots
import helper_fns

# return a dictionary holding eye limits from the 2021 and 2022 hurricane seasons!
# i decided to make a new file to simplify the analysis using the new tc datasets :)
def all_metadata():
    eyewall_time_limits = {}

    # part 1: get times from the 2021 eye pass data
    # look under "2023-02-27 get 2021 eyewall limits.ipynb" to see how I wrote some automatic code
    # to find these limits!
    eyewall_time_limits[ '2021'] = {}

    # date and time eye pass metadata
    dates = ['08-12-am', '08-12-pm', '0813',
            '08-16', '08-17', '08-18',
            '0819', '08-20', '0821',
            '08-27', '0829', '09-26',
            '09-27', '0929']
    passes = [
        [ (10.4787, 10.5821), (10.9949, 11.106), (12.3615, 12.5788)],
        [(22.1608, 22.3563), (23.3269, 23.4141)],
        [(9.6456, 9.8534), (10.9701, 11.194)],
        [(10.9583, 10.9938), (12.3583, 12.4216)],
        [(10.0851, 10.369), (11.4268, 11.5113), (12.868, 12.9924)],
        [(21.8553, 21.9103), (23.1903, 23.2798), (24.3673, 24.4395)],
        [(22.0095, 22.0995), (22.9457, 23.0707), (24.0498, 24.2193)],
        [(21.9232, 22.0104), (24.5149, 24.581)],
        [(22.5069, 22.5875), (23.7514, 23.823), (25.9095, 26.0095)],
        [(21.519, 21.669), (23.5469, 23.6558), (24.7547, 24.8097), (26.7216, 26.8372)],
        [(18.9966, 19.0349)],
        [(22.6209, 22.6826), (23.8448, 23.8937), (24.8804, 24.9521)],
        [(23.3316, 23.3971), (24.5296, 24.5724)],
        [(21.5788, 21.6793), (22.7988, 22.8877)]
        ]
    # add metadata to the dictionary!
    for datei, dateval in enumerate( dates):
        eyewall_time_limits['2021'][ dateval] = passes[ datei]

    # part 2: add 2022 eye pass data
    eyewall_time_limits[ '2022'] = {}

    dates = []
    passes = [

    ]





    return eyewall_time_limits
