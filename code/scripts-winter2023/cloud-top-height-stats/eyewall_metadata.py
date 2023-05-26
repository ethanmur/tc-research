import os
import numpy as np
import xarray as xr

os.chdir("/Users/etmu9498/research/code/scripts-winter2023")
import helper_fns_winter2023



#########
## IMPORTANT: should include eyewall case from 8/28 data!! there's a surprisingly
## clear eyewall pass that's been overlooked here!!!
#########



# return a dictionary holding eye limits from the 2021 and 2022 hurricane seasons!
# i decided to make a new file to simplify the analysis using the new tc datasets :)
def all_metadata( eye_limits='default'):
    eyewall_time_limits = {}

    # part 1: get times from the 2021 eye pass data
    # look under "2023-02-27 get 2021 eyewall limits.ipynb" to see how I wrote some automatic code
    # to find these limits!
    eyewall_time_limits[ '2021'] = {}
    eyewall_time_limits[ '2021']['eyewall_limits'] = {}
    eyewall_time_limits[ '2021']['intensity'] = {}
    eyewall_time_limits[ '2021']['category'] = {}
    eyewall_time_limits[ '2021']['names'] = {}
    # eyewall_time_limits[ '2021']['dates'] = {}


    # date and time eye pass metadata
    dates2021 = ['0811', '0812am', '0812pm', '0813',
            '0816', '0817', '0818',
            '0819', '0820', '0821',
            '0827', '0828', '0829',
            '0925', '0926', '0927', '0929']

    names2021 = [
        'fred', 'fred', 'fred', 'fred', 
        'grace', 'grace', 'grace', 
        'henri', 'henri', 'henri', 
        'ida', 'ida', 'ida', 
        'sam', 'sam', 'sam', 'sam', 
        ]

    # original eye pass limits! taken from old eyewall lims, matches previous results really well
    passes2021 = [
        [ ()], # no eyes for 8/11
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
        [ (21.165, 21.265)], # actually add an eye for 8/28!! one defintely exists, see the crl data...
        [(18.9966, 19.0349)],
        [()], # no eyes for 9/25
        [(22.6209, 22.6826), (23.8448, 23.8937), (24.8804, 24.9521)],
        [(23.3316, 23.3971), (24.5296, 24.5724)],
        [(21.5788, 21.6793), (22.7988, 22.8877)]
        ]
    # new eye pass limits. They remain the same for intense cases with well defined eyewalls,
    # but are at the -50 to 50 km limits for weaker cases. Inputed manually using scripts
    # found in notebooks-winter2023/"2022 eyewall determinations.ipynb"

    # 3/31/23 edit: sam 9/26 pass 0 has been shifted to the left
    passes2021_v2 = [
        [ ()], # no eyes for 8/11
        [ (10.329, 10.52), (11.049, 11.17), (12.55, 12.76)],
        [(22.092, 22.30), (23.348, 23.553)],
        [(10.057, 10.127), (10.98, 11.186)],
        [(10.89, 11.104), (12.17, 12.35)], # 8/16
        [(10.0851, 10.369), (11.315, 11.51), (12.79, 13.0)], # only pass 2, 3 changed
        [(21.8553, 21.9103), (23.1903, 23.2798), (24.3673, 24.4395)],
        [(22.0095, 22.0995), (22.9457, 23.0707), (24.0498, 24.2193)],
        [(21.775, 22.0), (24.415, 24.64)], # 8/20
        [(22.5069, 22.5875), (23.7514, 23.823), (25.9095, 26.0095)],
        [(21.519, 21.669), (23.5469, 23.6558), (24.7547, 24.8097), (26.7216, 26.8372)],
        [ (21.165, 21.265)], # actually add an eye for 8/28!! one defintely exists, see the crl data...
        [(18.9966, 19.0349)],
        [()], # no eyes for 9/25
        [(22.580, 22.6826), (23.8448, 23.8937), (24.8804, 24.9521)],
        [(23.3316, 23.3971), (24.5296, 24.5724)],
        [(21.5788, 21.6793), (22.7988, 22.8877)]
        ]

    # add 2021 intensity metadata! units: kt
    intensity = [
        0, 30, 32.5, 30,
        32.5, 45, 70,
        55, 60, 65,
        70, 85, 130,
        130, 135, 105, 115
        ]
    # add intensity and category metadata to the dictionary!
    for datei, dateval in enumerate( dates2021):
        eyewall_time_limits['2021']['intensity'][ dateval] = intensity[ datei]
    category = helper_fns_winter2023.find_category( intensity)
    for datei, dateval in enumerate( dates2021):
        eyewall_time_limits['2021']['category'][ dateval] = category[ datei]
    for datei, dateval in enumerate( names2021):
        eyewall_time_limits['2021']['names'][dateval] = names2021[datei]


    # part 2: add 2022 eye pass data
    eyewall_time_limits[ '2022'] = { }
    eyewall_time_limits[ '2022']['eyewall_limits'] = {}
    eyewall_time_limits[ '2022']['intensity'] = {}
    eyewall_time_limits[ '2022']['category'] = {}
    eyewall_time_limits[ '2022']['names'] = {}
    # eyewall_time_limits[ '2021']['dates'] = {}



    # working on creating 2022 metadata on this line
    # eyewall_time_limits[ '2022'] = {}
    dates2022 = [
         '0830', '0831', '0901', '0903', '0904', '0905', '0906', '0908',
         '0916', '0917', '0918', '0920', 
         '0924', '0925', '0926', '0927',
         '1007', '1008'
        ]
    names2022 = [
        'earl', 'earl', 'earl', 'earl', 'earl', 'earl', 'earl', 'earl', 
        'fiona', 'fiona', 'fiona', 'fiona', 
        'ian', 'ian', 'ian', 'ian',
        'julia', 'julia'
        ]

    # these passes use the original 2021 methods for finding cloud tops-
    # find the wind speed / cloud centers, then break things off at the 'eyewall'.
    # works well for intense cases, but likely misidentiifes eyewalls for weak
    # cases (especially after seeing radial distance plots!!).

    # try to develop a more quantitative approach for finding eye passes here.
    # only keep data 50 km away from the TC center for weak cases!
    # otherwise, use the same criteria for finding eyewalls
    passes2022 = [
          [ ( )], # 0830, na
          [ ()], # 8/31
          [ (9.13, 9.34), (10.42, 10.62), (12.6, 12.83)], # 9/01
          [ (9.28, 9.53), (10.17, 10.38), (11.35, 11.57), (11.64, 11.85)], # last case is too far -> (13.27, 13.52)], # 9/03
          [ (13.07, 13.3)], # 9/04
          [ (10.06, 10.32), (10.81, 11.02), (12.045, 12.193), (12.984, 13.201)], #09/05
          [ (11.18, 11.38), (12.29, 12.54), (13.365, 13.58)], # 0906, mid
          [ (10.44, 10.58), (11.79, 11.95), (12.93, 13.02)], # 0908, mid note!! crl and in situ times are off?!?
          [ (10.29, 10.514), (11.62, 11.83), ( 12.96, 13.196)], # 09/16
          [ (9.40, 9.60), (11.96, 12.197)],
          [ (9.97, 10.03), (11.2, 11.27), (12.38, 12.53), (13.27, 13.55), (14.29, 14.45)], # 0918, strong
          [ (10.37, 10.47), (11.43, 11.47), (12.77, 12.81), (13.96, 14.04)], #0920, strong
          [ (10.495, 10.81), (12.2, 12.4), (12.8, 13.) ], # 0924, na
          [ (10.67, 10.914), (11.884, 12.093)], # 09/25
          [ (10.15, 10.24), (11.42, 11.51), (11.60, 11.67), (12.69, 12.83)], # 0926, mid
          [ (10.38, 10.42), (11.7, 11.75), (12.58, 12.61)], # 9/27
          [ (22.4, 22.58), (23.45, 23.64), (24.65, 24.78)], # 10/07
          [ (22.1, 22.24), (23.16, 23.3), (24.9, 25.0)] # 1008, na
         ]
    # these passes were created using a new eye definition! Instead of looking at cloud structures,
    # the TC radial distance center is prioritized, with manual cutoffs being created around that.
    # 3/14 update: not really using this case...
    passes2022_v2 = [
          [ ( )], # 0830, na
          [ (9.80, 10.18), (10.90, 11.15), (12.15, 12.3), (12.95, 13.3)], # 0905, weak
          [ (11.22, 11.40), (12.30, 12.56), (13.41, 13.60)], # 0906, mid
          [ (10.45, 10.57), (11.79, 11.95), (12.93, 13.02)], # 0908, mid note!! crl and in situ times are off?!?
          [ (10.30, 10.47), (11.58, 11.74), (13.08, 13.19)], # 0916, weak
          [ (9.56, 9.78),  (11.92, 12.17)], # 0917, weak
          [ (9.93, 10.05), (11.17, 11.28), (12.38, 12.49), (13.27, 13.57), (14.30, 14.48)], # 0918, strong
          [ (10.37, 10.44), (11.43, 11.49), (12.78, 12.85), (13.95, 14.03)], #0920, strong
          [ ()], # 0924, na... maybe add some??
          [ (10.72, 10.9), (12.02, 12.10)], # 0925, weak
          [ (10.18, 10.25), (11.42, 11.51), (11.60, 11.67), (12.69, 12.83)], # 0926, mid
          [ ()] # 1008, na
         ]

    # set 2021 and 2022 eyewalls depending on user input!
    if eye_limits == 'passes':
        for datei, dateval in enumerate( dates2021):
            eyewall_time_limits['2021']['eyewall_limits'][ dateval] = passes2021[ datei]
        for datei, dateval in enumerate( dates2022):
            eyewall_time_limits['2022']['eyewall_limits'][ dateval] = passes2022[ datei]
    elif eye_limits == 'passes2':
        for datei, dateval in enumerate( dates2022):
            eyewall_time_limits['2022']['eyewall_limits'][ dateval] = passes2022_v2[ datei]

    # this is the default case! use the updated eyewall limits
    elif eye_limits == 'default' or eye_limits == 'passes 50km':
        for datei, dateval in enumerate( dates2021):
            eyewall_time_limits['2021']['eyewall_limits'][ dateval] = passes2021_v2[ datei]
        for datei, dateval in enumerate( dates2022):
            eyewall_time_limits['2022']['eyewall_limits'][ dateval] = passes2022[ datei]

    # units: kt
    # first tc is irrelevant test data
    # data manually pulled from the noaa advisory archive https://www.nhc.noaa.gov/archive/2022/
    intensity = [
        25, 25, 25, 40, 45, 55, 55, 85, # earl
        50, 50, 55, 100, # fiona 
        35, 40, 70, 110, # ian
        45, 65 # julia
        ]
    # add intensity, name, and category metadata to the dictionary!
    for datei, dateval in enumerate( dates2022):
        eyewall_time_limits['2022']['intensity'][ dateval] = intensity[ datei]
    category = helper_fns_winter2023.find_category( intensity)
    for datei, dateval in enumerate( dates2022):
        eyewall_time_limits['2022']['category'][ dateval] = category[ datei]
    for datei, dateval in enumerate( names2022):
        eyewall_time_limits['2022']['names'][ dateval] = names2022[ datei]

    return eyewall_time_limits
