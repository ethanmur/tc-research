########################
##   run_radQC        ##
##   VERSION: 1.3.8   ##
########################

# To run this program and display all output (stdout and stderr) on the console only:
#
#       python run_radQC.py input_radQC.yml
#
# To run in the background and write all output to an output file (noting that you can
#   provide any path/filename you want for the output file):
#       
#       python -u run_radQC.py input_radQC.yml &> /path/and/filename/desired/for/crqc.out &
#
# To run while simultaneously writing to the console AND some output file:
#
#       python -u run_radQC.py input_radQC.yml 2>&1 | tee /path/and/filename/desired/for/crqc.out



# Below are output variable names, including those output when output of intermediate
# versions of edited variables is enabled.
#
# Intermediate variable versions use the output reflectivity and velocity names as a base.
# For example, if the output reflectivity (`refOut`) and velocity (`velOut`) variables are 
# named 'DD' and 'VD', respectively, intermediate variables will all begin with 'DD_' 
# and 'VD_'.
#
#
# 
# Output when intermediate version output is enabled:
# <refOut>_0A:    Reflectivity prior to running spoke removal routine
# <refOut>_0spkD: Reflectivity difference from previous gate (in azimuth at constant range) when a potential spoke is detected
# <refOut>_0sDd:  Reflectivity difference between <refOut>_0A and <refOut>_0B (with masked/missing values filled with zeros)
# <refOut>_0B:    Reflectivity after running spoke removal routine
# ----------
# <refOut>_1A:    Reflectivity prior to running surface removal routine
# Z_cSum:         Cumulative sum of Z (derived from <refOut>_1A) (normalized by 1x10^6) along each ray
# <refOut>_1g:    Reflectivity following geometric surface removal
# <refOut>_1gDd:  Reflectivity difference between <refOut>_1A and <refOut>_1g (with masked/missing values filled with zeros)
# <refOut>_1ag:   Mean + Std. dev. of aggregate reflectivity field
# <refOut>_1mDd:  Reflectivity difference between <refOut>_1g and <refOut>_1B (with masked/missing values filled with zeros)
# <refOut>_1B:    Reflectivity after running surface removal routine
# ----------
# <refOut>_2A:    Reflectivity prior to running 2D gaussian filter
# <refOut>_2sm:   Smoothed reflectivity field.
# <refOut>_2fDd:  Reflectivity difference between <refOut>_2A and <refOut>_2B (with masked/missing values filled with zeros)
# <refOut>_2B:    Reflectivity after running 2D gaussian filter
# ----------
# <velOut>_0A:    Velocity prior to running spoke removal routine
# <velOut>_0B:    Velocity after running spoke removal routine
# ----------
# <velOut>_1A:    Velocity prior to running surface removal routine
# <velOut>_1g:    Velocity following geometric surface removal
# <velOut>_1B:    Velocity after running surface removal routine
# ----------
# <velOut>_2A:    Velocity prior to running 2D gaussian filter
# <velOut>_2B:    Velocity after running 2D gaussian filter
# ----------
# <velOut>_3A:    Velocity prior to running clear air thresholding using an environmental wind profile
# <velOut>_3B:    Velocity after running clear air thresholding using an environmental wind profile

# Convert to CF radial first
import warnings
warnings.filterwarnings("ignore", message='Velocities outside of the Nyquist interval found in sweep', append=True)
warnings.filterwarnings("ignore", message='All-NaN slice encountered', append=True)
warnings.filterwarnings("ignore", message='invalid value encountered in sqrt', append=True)
warnings.filterwarnings("ignore", message='overflow encountered in power', append=True)
warnings.filterwarnings("ignore", message='invalid value encountered in less', append=True)

import pyart
import sys
import os
import gc
from os import environ as _environ
from subprocess import call, run
from glob import glob
import numpy as np
from datetime import datetime as dt
import datetime
from shutil import copy2
import pandas as pd
import yaml

from radQC_funcs import *




inputFile = sys.argv[1]


## ** New input method - DS (10/28/19) ** ##
with open(inputFile, 'r') as pydInPrms:
    prmsIn = yaml.load(pydInPrms,Loader=yaml.FullLoader)

for key,val in prmsIn.items():
        exec(key + '=val')


fmtList = ['dorade','foray','nc']
if frmat == 'None':
    pass
elif frmat not in fmtList:
    print('Please specify one of the following formats...')
    print(fmtList)
    quit()

# Set the CRQC_QUIET environment variable to suppress authorship/citation info on run
#   e.g., if you're using the bash shell:
#       export CRQC_QUIET=1
#   You may also wish to suppress the citation block for Py-ART:
#       export PYART_QUIET=1
if 'CRQC_QUIET' not in _environ:
    print("""
    ##################################################################################

    Authors: Daniel Betten and Addison Alford

    This script dealiases radial velocity, differential phase, specific differential
    phase, reflectivity, and differential reflectivity using PyART as the primary tool
    (Helmus and Collins 2016, JORS). As of 23 February 2017, differential phase and
    specific differential phase are smoothed and recalculated using a least squares
    fit to a 3 km radial window. Attenuation is corrected using a self-consistent Z-Phi
    method (Gu et al. 2011, JAMC). The "strong rotation" option available in the input
    file should not be used as of yet. This will be updated with time. All questions
    can be directed to addisonalford@ou.edu.

    ##################################################################################

    ** Ground target removal and filtering**
    Author: Dan Stechman

    This addition to the script provides methods for identifying and masking ground
    targets, specifically for use with data collected by helically-scanning radars
    such as the NOAA Tail Doppler Radar (TDR). Additional methods for reducing noise
    and other artifacts are also provided. Questions pertaining to this portion of the
    software suite can be directed to daniel.stechman@noaa.gov.

    ##################################################################################

    """)

    if (stag) or (dual):
        print("""
        ##################################################################################
        The staggered-PRT & dual-PRF processing error correction was written by
        Addison Alford, and the result of collaboration between Dr. Michael Biggerstaff
        Curtis Riganti, and Gordon Carrie. The software is still a work in progress,
        and any changes or suggestions for improvement are welcomed. Please send
        any suggestions to addisonalford@ou.edu.
    
        The dual-PRF processing is a modified from Joe and May (2003, JTECH).
    
        ##################################################################################
    
        """)

print('\nCurrent date/time: {:%Y/%m/%d %H:%M:%S}\n'.format(dt.now()))

print('********************************************************************')
print('*************** User-defined parameters for this run ***************')
print('********************************************************************\n')

print('Input CfRadial directory = {}'.format(radar_folder))
print('Platform is airborne = {}'.format(air))
print('Input velocity field name = {}'.format(velIn))
print('Input reflectivity field name = {}'.format(refIn))
print('Output velocity field name = {}'.format(velOut))
print('Output reflectivity field name = {}'.format(refOut))

print('\nSave intermediate versions of output variables = {}'.format(svExtra))
print('Output format = {}'.format(frmat))

print('\nRemove "spokes" prior to aggregation/editing = {}'.format(rmvSpokes))
if rmvSpokes:
    print(' > DBZ diff threshold to flag a gate as a spoke candidate = {} dBZ'.format(spokeDiff))
    print(' > Range beyond which to mask ALL data (0 disables) = {} m'.format(maskRingOut))
    print(' > Range within which ALL data will be masked (0 disables) = {} m'.format(maskRingIn))
    print(' > Apply dual-PRT reflectivity offset correction = {}'.format(dualPRTofst))

print('\nGenerate PDD aggregate stats for airborne ground echo removal = {}'.format(agg))
if agg:
    if vols[0][1] == 0 or vols[0][2] == 0:
        if vols[0][3] == 'L':
            posStr = 'left'
        elif vols[0][3] == 'R':
            posStr = 'right'
        print(' > Volume set by user to include all sweep files in radar_folder.')
        print('     Storm is to the {} of the aircraft'.format(posStr))
    else:
        print(' > Volumes defined in input file:')
        print('   Volume\tStart\t\t\tEnd\t\t\tStorm Pos.')
        print('   -----------------------------------------------------------------------')
        for vol in vols:
            print('   {}\t\t{:%Y/%m/%d %H:%M:%S}\t{:%Y/%m/%d %H:%M:%S}\t{}'.format(vol[0],dt.strptime(vol[1],'%Y%m%d_%H%M%S'),
                                                                                   dt.strptime(vol[2],'%Y%m%d_%H%M%S'),vol[3]))

print('\nGenerate simulated radial velocity from input wind profile = {}'.format(makeSimVel))
if makeSimVel:
    print(' > Input wind profile = {}'.format(soundF))
    print(' > Run QC and smoothing on wind profile = {}'.format(qcSnd))
    print(' > Name to assign to the output simulated radial velocity field = {}'.format(simVelID))

print('\nDealias radial velocity = {}'.format(deal))
if deal:
    print(' > Dealias with the 4DD algorithm using an input environmental wind profile = {}'.format(deal4d_snde))
    print(' > Dealias with the 4DD algorithm using a previously edited volume = {}'.format(deal4d_vol))
    if deal4d_vol:
        print('   > Previously edited volume = {}'.format(editVol))
    print(' > Perform 2nd dealiasing pass (only if deal4d_* options are disabled) = {}'.format(dealias2))
    if dealias2:
        print('   > Strong rotation present = {}'.format(strong_rot))
    print(' > Exclude a given range of azimuths from dealiasing = {}'.format(azRstrctDeal))
    if azRstrctDeal:
        print('   > Min azimuth to exclude = {}°'.format(minAzExcld))
        print('   > Max azimuth to exclude = {}°'.format(maxAzExcld))
    print(' > SNR field name to use when dealiasing = {}'.format(snrn))
    print('   > SNR field is actually SNR = {}'.format(SNR))
    print(' > Trip flag field name = {}'.format(flg))
    print(' > SNR field threshold to mask below = {}'.format(snrThres))
    print(' > Velocity standard dev. to mask above = {}'.format(stdThresh))
    print('   > Window length over which to calculate the std. dev. = {}'.format(windowLen))
    print('   > Max reflectivity to apply std. dev. filter = {} dBZ'.format(refThresh))
    print(' > Number of gates used when identifying speckles (0 disables) = {}'.format(despk))

print('\nCorrect dual-PRF processor mistakes = {}'.format(dual))
if dual:
    print(' > Use new version of dual-PRF code = {}'.format(useNew_dPRF))

print('\nCorrect staggered-PRT processor mistakes = {}'.format(stag))

print('\nKDP calculation using spline (fast) = {}'.format(kdpCalcFast))
print('KDP calculation using least squares (slow) = {}'.format(kdpCalcSlow))
if kdpCalcFast or kdpCalcSlow:
    print(' > PHIDP field name to use in KDP calculation = {}'.format(phin))
    print(' > KDP field name to use in KDP calculation = {}'.format(kdpn))
    print(' > RhoHV field name to use in KDP calculation = {}'.format(rhon))

print('\nPerform attenuation correction = {}'.format(attenCorr))
if attenCorr:
    print(' > ZDR field name to use in attentuation correction = {}'.format(zdrn))

print('\nMask ground returns from airborne radar data = {}'.format(sfcRmv))
if sfcRmv:
    print(' > Radar altitude adjustment (only for geometric surface masking) = {} m'.format(altAdj))
    print(' > Aggregate DBZ mean + std. dev. above which to consider masking gates = {} dBZ'.format(aggMagThrsh))
    print(' > Max absolute value of velocity expected for ground targets = {} m/s'.format(sfcVel))
    print(' > Degrees (rotation angle) above radar-relative horizon to include in surface masking (storm side) = {}°'.format(mskAzBuf_strm))
    print(' > Degrees (rotation angle) above radar-relative horizon to include in surface masking (clear side) = {}°'.format(mskAzBuf_clr))
    print(' > Threshold cumulative sum of Z along radial below which to apply surface masking = {} (× 1×10^6 mm^6 m^-3)'.format(zCsumThrsh))

print('\nReduce noise using 2D Gaussian filter = {}'.format(g2dFilt))
if g2dFilt:
    print(' > Run filtering prior to dealiasing (True) or as final editing step (False) = {}'.format(preDealThrsh))
    print(' > Variable to use for thresholding = {}'.format(thrshVarID))
    print(' > Value to fill missing/NaN gates with in thresholding field prior to smoothing = {}'.format(thrshFillVal))
    print(' > Std. dev. value in X used in filter = {}'.format(filtStdX))
    print(' > Std. dev. value in Y used in filter = {}'.format(filtStdY))
    print(' > Value of smoothed thresholding field to threshold below = {}'.format(smthThrsh))
    print(' > Turn off filtering for short pulse gates = {}'.format(retainShrtPulse))

print('\nMask clear air velocities using simulated radial velocity field = {}'.format(clrVelChk))
if clrVelChk:
    print(' > Max value of reflectivity to be considered clear air = {} dBZ'.format(refThrsh))
    print(' > Absolute difference allowed between the observed and simulated radial velocity before masking = {} m/s'.format(velDifThrsh))

if sfcRmv and not air:
    print('\n\t>>> WARNING <<< Surface removal was designed for airborne radar. Airborne' 
          '\n\toption was disabled in input file. Exiting...')
    sys.exit()


print('********************************************************************\n')

# Save a copy of the input YAML file to the output directory, prepended by
#   the current date/time. This may be useful for determining later on exactly 
#   what options were enabled and what parameters were used.
# if svPrms:
#     ymlF_log = copy2(inputFile,'{}/_{:%Y%m%d-%H%M%S}_inputPrms.yml'.format(radar_folder,dt.now()))


wholeTime0 = dt.now()


if rmvSpokes:
    fileList = np.asarray(sorted(glob(radar_folder+'/cfrad.*')))

    print('- Removing spokes from radar fields...')
    
    for i,radF in enumerate(fileList):
        print('\tSweep {} of {}...'.format(i+1,len(fileList)))
        radar = pyart.io.read_cfradial(radF)
        
        if velOut not in radar.fields.keys():
            radar.add_field_like(velIn,velOut,radar.fields[velIn]['data'].copy())
        elif velOut == velIn:
            print('\tvelOut and velIn are the same. Adjust `velOut` in input YAML file to avoid overwriting input data. Exiting...')
            sys.exit()
        if refOut not in radar.fields.keys():
            radar.add_field_like(refIn,refOut,radar.fields[refIn]['data'].copy())
        elif refOut == refIn:
            print('\trefOut and refIn are the same. Adjust `refOut` in input YAML file to avoid overwriting input data. Exiting...')
            sys.exit()
        
        if svExtra:
            radar.add_field_like(refOut,refOut + '_0A',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_0A',radar.fields[velOut]['data'].copy(),replace_existing=True)
        
        radar = rmvSpoke(radar,refOut,velOut,air,spokeDiff=spokeDiff,maskRingOut=maskRingOut,
                         maskRingIn=maskRingIn,dualPRTofst=dualPRTofst,svExtra=svExtra)
        
        if svExtra:
            radar.add_field_like(refOut,refOut + '_0B',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_0B',radar.fields[velOut]['data'].copy(),replace_existing=True)
    
        radF_new = radF[:-3]+'_preCorr.nc'
        pyart.io.write_cfradial(radF_new,radar)


if agg or sfcRmv:
    if rmvSpokes or any('pre' in os.path.basename(rF) for rF in sorted(glob(radar_folder+'/cfrad.*'))):
        initFileList = sorted(glob(radar_folder+'/*_preCorr.nc'))
        aggFileList = np.asarray(initFileList)
    else:
        initFileList = sorted(glob(radar_folder+'/cfrad.*'))
        aggFileList = np.asarray([fn for fn in initFileList if ('_corr' not in os.path.basename(fn))])

    filesDTstr = [os.path.split(fN)[1][6:21] for fN in aggFileList]
    filesDT = np.asarray([dt.strptime(aDT,'%Y%m%d_%H%M%S') for aDT in filesDTstr])

    
    # Pull apart the nested list with volume start/end times and storm positions
    volNum = [x[0] for x in vols] # Not currently used - included to make visual ID of volumes easier in input file
    volStrtDT = [x[1] for x in vols]
    volEndDT = [x[2] for x in vols]
    volStrmPos = [x[3] for x in vols]


    if volStrtDT == [0]:
        volStrtT = [filesDT[0]]
    else:
        volStrtT = [dt.strptime(sT,'%Y%m%d_%H%M%S') for sT in volStrtDT]
        
    if volEndDT == [0]:
        volEndT = [filesDT[-1]]
    else:
        volEndT = [dt.strptime(eT,'%Y%m%d_%H%M%S') for eT in volEndDT]
        

if agg:
    outFnAgg = []
    for i,(iVs,iVe,iVp) in enumerate(zip(volStrtT,volEndT,volStrmPos)):
        print('- Creating mask for volume spanning {:%Y/%m/%d %H:%M:%S} to {:%Y/%m/%d %H:%M:%S}...'.format(iVs,iVe))
        
        volMatchIx = np.squeeze(np.where(np.logical_and(filesDT >= iVs, filesDT <= iVe)))
        filesAgg = aggFileList[volMatchIx]
    
        if len(filesAgg) > 0:
            radAgg = genAgg(filesAgg,iVp,velOut,refOut,velIn,refIn)
        else:
            print('\tNo cfradial files were found in radar_folder within the specified volume\n' 
                  '\tstart and end times. Exiting...')
            sys.exit()
        
        swpAngl = radAgg.fixed_angle['data'][0]
    
        if swpAngl < 0:
            swpPnt = 'AFT'
        else:
            swpPnt = 'FORE'

        print('\t*** Saving aggregate fields...')
        #tmpRadar = pyart.io.read_cfradial(aggFileList[0])
        #strtFStr = pyart.graph.common.generate_radar_time_begin(tmpRadar).strftime("%Y%m%d_%H%M%S")
        #tmpRadar = pyart.io.read_cfradial(aggFileList[-1])
        #endFStr = pyart.graph.common.generate_radar_time_begin(tmpRadar).strftime("%Y%m%d_%H%M%S")
        strtFStr = iVs.strftime("%Y%m%d_%H%M%S")
        endFStr = iVe.strftime("%Y%m%d_%H%M%S")
        
        outFnAgg.append('{}/aggCfrad.{}_{}_to_{}.nc'.format(radar_folder,swpPnt,strtFStr,endFStr))
        pyart.io.write_cfradial(outFnAgg[i],radAgg)


if not agg and sfcRmv:
    outFnAgg = sorted(glob(radar_folder+'/aggCfrad.*'))

    if len(outFnAgg) != len(volStrtT):
        print('Number of aggCfrad files does not match number of volume start/end times. Exiting...')
        sys.exit()

if sfcRmv:        
    if len(outFnAgg) == 1:
        applyAggStrtT = [filesDT[0]]
        applyAggEndT = [filesDT[-1]]
    
    else:
        applyAggStrtT = []
        applyAggEndT = []
        for iV in range(len(volStrtT)):
            if iV == 0:
                applyAggStrtT.append(filesDT[0])
                applyAggEndT.append(volStrtT[iV+1] - datetime.timedelta(seconds=1))
            elif iV != len(volStrtT) - 1:
                applyAggStrtT.append(volStrtT[iV])
                applyAggEndT.append(volStrtT[iV+1] - datetime.timedelta(seconds=1))
            else:
                applyAggStrtT.append(volStrtT[iV])
                applyAggEndT.append(filesDT[-1])
                


if makeSimVel or (deal and deal4d_snde) or clrVelChk:
    import metpy.calc as mpcalc
    from metpy.units import units
    
    if not makeSimVel:
        print('!! makeSimVel was set to False, but deal4d_snde and clrVelChk require the associated output. Overriding and continuing...')
    print('- Generating environmental wind profile from sounding...')
    print('\tReading environmental sounding from {}'.format(soundF))
    
    sndDat = pd.read_csv(soundF,header=[12],delim_whitespace=True,skiprows=[13,14])
    sndDat.drop_duplicates(subset='Alt',inplace=True)
    if qcSnd:
        print('\tRunning QC and smoothing on profile')
        sndDat['Alt'] = sndDat['Alt'].replace(99999.0,np.nan)
        sndDat[['Ucmp','Vcmp']] = sndDat[['Ucmp','Vcmp']].replace(9999.0,np.nan)
        sndDat[['Qu','Qv']] = sndDat[['Qu','Qv']].replace(99.0,np.nan)
        
        dropIx = sndDat[(sndDat['Qu'] != 1.0) | (sndDat['Qv'] != 1.0)].index
        sndDat.drop(dropIx,inplace=True)

        # Drop any rows with all NaN values for winds
        sndDat_clean = sndDat.dropna(subset=('Ucmp','Vcmp'), how='all').reset_index(drop=True)
        u1 = sndDat_clean['Ucmp'].values * units.meter_per_second
        v1 = sndDat_clean['Vcmp'].values * units.meter_per_second
        
        u = mpcalc.smooth_gaussian(u1,8)
        v = mpcalc.smooth_gaussian(v1,8)
        alt_msl = sndDat_clean['Alt'].values
    
    else:
        u = sndDat['Ucmp'].values
        v = sndDat['Vcmp'].values
        alt_msl = sndDat['Alt'].values
    
    if len(alt_msl) > 999:
        print('\tToo many observations in sounding ({}) - max is 999. Resampling...'.format(len(alt_msl)))
        # Set spacing interval--Every 20 meters from 0 to 19000 meters
        intvl = np.arange(0, 19000, 20) * units.meter
        
        # Get indexes of values closest to defined interval
        rsmplIx = mpcalc.resample_nn_1d(alt_msl * units.meter, intvl)
        
        u = u[rsmplIx]
        v = v[rsmplIx]
        alt_msl = alt_msl[rsmplIx]
        
    
    envWndProf = pyart.core.HorizontalWindProfile.from_u_and_v(alt_msl,u,v)


filePrev = editVol
errorList = []

if rmvSpokes or any('pre' in os.path.basename(rF) for rF in sorted(glob(radar_folder+'/cfrad.*'))):
    radFileList = sorted(glob(radar_folder+'/*_preCorr.nc'))
else:
    radFileList = sorted(glob(radar_folder+'/cfrad.*'))


radFiles = [fn for fn in radFileList if ('_corr' not in os.path.basename(fn)) and ('RHI' not in os.path.basename(fn))]


for i,radar_file in enumerate(radFiles):

    print('\n********************************************************************')
    print('** Working on file {} of {}...'.format(i+1,len(radFiles)))    
    
    time0 = dt.now()
    
    if radar_file == filePrev:
        if deal4d_vol:
            call(['RadxConvert -%s -outdir %s -disag -f %s' % (frmat,radar_file[:radar_file.rindex('cf')-1],radar_file)],shell=True)
        continue
        
    print('Editing {}'.format(radar_file))
    
    radar = pyart.io.read_cfradial(radar_file)
    
    
    elev = round(radar.get_elevation(0).mean(), 2)
    if elev > 100.:
        continue
    
    nyq = round(radar.get_nyquist_vel(0), 3)
    radar.instrument_parameters['nyquist_velocity']['data'][:] = nyq
    
    
    
    # Initialize output velocity/reflectivity variables if they don't already exist
    if velOut not in radar.fields.keys():
        radar.add_field_like(velIn,velOut,radar.fields[velIn]['data'].copy())
    elif velOut == velIn:
        print('velOut and velIn are the same. Adjust `velOut` in input YAML file to avoid overwriting input data. Exiting...')
        sys.exit()
    if refOut not in radar.fields.keys():
        radar.add_field_like(refIn,refOut,radar.fields[refIn]['data'].copy())
    elif refOut == refIn:
        print('refOut and refIn are the same. Adjust `refOut` in input YAML file to avoid overwriting input data. Exiting...')
        sys.exit()
    
    
    if makeSimVel or (deal and deal4d_snde) or clrVelChk:
        print('- Generating simulated radial velocity field...')
        simVel = pyart.util.simulated_vel_from_profile(radar, envWndProf, 
                                                       sim_vel_field=simVelID,airTail=air)
        radar.add_field(simVelID,simVel)
    
    # Perform an optional thresholding operation prior to attempting dealiasing
    # Doing this tends to significantly speed up and improve dealiasing
    if g2dFilt and preDealThrsh:
        print('- Running pre-dealiasing noise filtering...')
    
        if svExtra:
            radar.add_field_like(refOut,refOut + '_2A',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_2A',radar.fields[velOut]['data'].copy(),replace_existing=True)
        
        radar = radFilt(radar,thrshVarID,refOut,velOut,thrshFillVal,
                        filtStdX=filtStdX,filtStdY=filtStdY,
                        smthThrsh=smthThrsh,retainShrtPulse=retainShrtPulse,
                        svExtra=svExtra)
                        
        if svExtra:
            radar.add_field_like(refOut,refOut + '_2B',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_2B',radar.fields[velOut]['data'].copy(),replace_existing=True)
        
    
    if deal:
        if deal4d_snde and deal4d_vol:
            print('When using the 4DD dealiasing method, currently only one of deal4d_snde and deal4d_vol can be true. Exiting...')
            sys.exit()
    
        if deal4d_snde:
            print('- Dealiasing radial velocity using the 4DD method (using input environmental sounding)...')
        elif deal4d_vol:
            print('- Dealiasing radial velocity using the 4DD method (using previous edited volume)...')
        else:
            print('- Dealiasing radial velocity using the region-based method...')

        if snrn not in radar.fields.keys():
            print('\t{} not found in radar file. Skipping...'.format(snrn))
            continue
    
        if radar.azimuth['data'].shape[0] > 358: 
            wrap = True
        else: 
            wrap = False
        #print('\tVolume is a 360 PPI = {}'.format(wrap))

        if stag or dual: 
            print('\tExtended Nyquist is {:.3f}'.format(nyq))
        else: 
            print('\tNyquist is {:.3f}'.format(nyq))
        radar.instrument_parameters['nyquist_velocity']['data'][:] = np.round(nyq,3)

        snr = radar.fields[snrn]['data'].copy()
        vel = radar.fields[velOut]['data'].copy()
        ref = radar.fields[refOut]['data'].copy()

        #Check to see if velocity fields exist.
        if flg != 'None': 
            trip = radar.fields[flg]['data'].copy()

        if 'VEL1' not in radar.fields.keys():
            radar.add_field('VEL1',radar.fields[velOut].copy())
        
        if SNR:
            ref = np.ma.masked_array(ref, snr == 0)
        
        if flg != 'None':
            vel = np.ma.masked_where(np.ma.masked_outside(\
                radar.fields[flg]['data'],0.,1.5).mask == True,vel)
            radar.fields[refOut]['data'] = np.ma.masked_where(\
                np.ma.masked_outside(radar.fields[flg]['data'],0.,1.5).mask == True,ref)

    
        #First threshold velocity field on SNR field below value of snrThres
        try:
            radar.fields[velOut]['data'] = std_filter(vel,snr,ref,refThresh = refThresh,
                                                    mthresh=snrThres,sthresh=stdThresh,
                                                    window=windowLen,snr=SNR)
            # (DS - 1/28/20) - Commented out as it's not clear why this version was used for only certain cases...
            #for sweep_slice in radar.iter_slice():
            #    radar.fields[velOut]['data'][sweep_slice] = std_filter(radar.fields[velOut]['data'][sweep_slice],
            #                                                         snr[sweep_slice],ref[sweep_slice],
            #                                                         mthresh=snrThres,sthresh=stdThresh,
            #                                                         window=windowLen,snr=SNR)
        except:
            print('\tSomething went wrong while thresholding on SNR.')
            traceback.print_exc()
            sys.exit()
            #continue
    
        #Despeckle the velocity field.
        despeckle(radar,vel_field=velOut, ngates=despk)
    
        if azRstrctDeal:
            V_preDea1 = radar.fields[velOut]['data'].copy()
            az_raw = radar.rotation['data'].data
            roll = radar.roll['data'].data
            azmth = az_raw + roll
            azmth[azmth < 0] += 360
            azmth[azmth > 360] -= 360
            azExcldIx = np.where(np.logical_and(azmth >= minAzExcld, azmth <= maxAzExcld))[0]
    
        #Get the mask and dealias the velocity field. To get a better representation
        #   of the nyquist, bump up the interval_splits. Default is 3. More equals
        #   separating the nyquist interval into more splits. This creates more 
        #   regions. The more regions, the longer the program takes. However,
        #   data should be better in the end.
        Vmask = radar.fields[velOut]['data'].mask.copy()
    
        radar.fields[velOut]['data'] = var_fill(radar.fields[velOut]['data'].copy())
    
    
        if deal4d_snde:
            radar.fields['VEL1'] = pyart.correct.dealias_fourdd(radar, sonde_profile=envWndProf, 
                                                    vel_field=velOut, corr_vel_field=velOut)
    
        elif deal4d_vol:
            if rhon not in radar.fields.keys():
                print('\t{} not found in radar file and is required for 4DD (volume) dealiasing. Exiting...'.format(rhon))
                sys.exit()
            gateFilt = pyart.filters.moment_based_gate_filter(radar,ncp_field=None,
                                                              rhv_field = rhon,refl_field=refOut,
                                                              min_ncp=0,min_rhv=0.02,
                                                              min_refl=-5,max_refl=100)
            
            radar.fields['VEL1']['data'] = np.ma.masked_where(gateFilt.gate_excluded == True,
                                                              radar.fields['VEL1']['data'])
        
            radar.fields[velOut]=pyart.correct.dealias_fourdd(radar, last_radar=radar_filePrev, 
                                                            sonde_profile=None, gatefilter=gateFilt,
                                                            filt=0, rsl_badval=131072.0, 
                                                            keep_original=False, set_limits=False,
                                                            vel_field=velOut, corr_vel_field=velOut, 
                                                            last_vel_field=velOut,debug=False, 
                                                            max_shear=0.05, sign=-1,proximity=5,
                                                            mingood=5,stdthresh = 0.5,thresh = 0.7,
                                                            compthresh = 0.5,compthresh2 = 0.6)
        else:
            radar.fields['VEL1'] = pyart.correct.dealias_region_based(radar,vel_field=velOut,
                                                    interval_splits=3,rays_wrap_around=wrap)
                                                
        if not deal4d_vol:
            radar.fields['VEL1']['data'] = np.ma.masked_array(radar.fields['VEL1']['data'],
                                                              Vmask==True)
            radar.fields[velOut]['data'] = np.ma.masked_array(radar.fields[velOut]['data'],
                                                            Vmask==True)
        
            despeckle(radar,vel_field='VEL1', ngates=despk)
        
        
        
            dea, deaF = varsigndiff(radar.fields[velOut]['data'].copy()[:,1:],
                                    radar.fields['VEL1']['data'].copy()[:,1:])
                                
            #Do a second pass on the corrected velocity field.
            if dealias2 and deal4d_snde:
                print('\tA second run of 4DD dealiasing is generally not required and will thus be skipped...')
            if dealias2 and not deal4d_snde:
                print('\tBeginning second pass of region-based dealiasing...')
                if azRstrctDeal:
                    radar.fields['VEL1']['data'][azExcldIx,:] = V_preDea1[azExcldIx,:]
                    V_preDea2 = radar.fields['VEL1']['data'].copy()
                placeholder = radar.fields['VEL1']['data'].copy()
                var1_fill = var_fill(placeholder)
    
                if 'VEL2' not in radar.fields.keys():
                    radar.add_field_like(velOut,'VEL2',radar.fields[velOut]['data'].copy())
    
                varg2 = dealias_corr(var1_fill[:,1:].copy(),deaF,nyq)
                radar.fields['VEL2']['data'][:,1:] = np.ma.masked_array(varg2,Vmask[:,1:] ==True)
                radar.fields['VEL3'] = pyart.correct.dealias_region_based(radar,vel_field='VEL2',
                                                                          interval_splits=6,
                                                                          rays_wrap_around=wrap)
    
    
                radar.fields[velOut]['data'] = radar.fields['VEL3']['data'].copy()
        

                vel = radar.fields[velOut]['data'].copy()
                radar.fields[velOut]['data'] = std_filter(vel,snr,ref,refThresh,
                                                        mthresh=snrThres,sthresh=stdThresh,
                                                        window=windowLen,snr=SNR)
                                                    
                #If there is strong rotation, do yet another pass.
                if strong_rot:
                    if 'VEL5' not in radar.fields.keys():
                        radar.add_field_like(velOut,'VEL5',radar.fields[velOut]['data'].copy())
                    dea2, deaF2 = varsigndiff(radar.fields[velOut]['data'].copy()[:,1:],
                                              radar.fields['VEL2']['data'][:,1:])
            
                    varg3 = var_rot_tel(radar.fields['VEL2']['data'][:,1:].copy(),dea2,deaF2,nyq)
            
                    radar.fields['VEL5']['data'][:,1:] = np.ma.masked_array(varg3,Vmask[:,1:] ==True)
            
                    despeckle(radar,vel_field='VEL5', ngates=despk)
            
                    radar.fields[velOut]['data'] =radar.fields['VEL5']['data']
            
                if azRstrctDeal:
                    radar.fields[velOut]['data'][azExcldIx,:] = V_preDea2[azExcldIx,:]
    
            elif not dealias2:
                radar.fields[velOut]['data'] = radar.fields['VEL1']['data'].copy()
        
                if azRstrctDeal:
                    radar.fields[velOut]['data'][azExcldIx,:] = V_preDea1[azExcldIx,:]

        
                vel = radar.fields[velOut]['data'].copy()
                radar.fields[velOut]['data'] = std_filter(vel,snr,ref,refThresh,
                                                        mthresh=snrThres,sthresh=stdThresh,
                                                        window=windowLen,snr=SNR)
                                                    
        #Do a final despeckle on all fields.
        despeckle(radar,vel_field=velOut,ngates=despk)
    
        gc.collect()
    
    
    #Calculation of 'folding Nyquist' for dual-PRF corrections.
    #vLow refers to the Nyquist of the low ratio (high PRF).
    if stag or dual:
        if stag:
            print('- Correcting staggered-prt errors...')
            print('\tExtended Nyquist is {:.3f}'.format(nyq))
        if dual:
            print('- Correcting dual-PRF errors...')
            print('\tExtended Nyquist is {:.3f}'.format(nyq))
            
        freq = radar.instrument_parameters['frequency']['data'].copy()[0]
        c = 2.9989e8
        wLen = c/freq
        
        if 'DOW' in radar_file:
            vLow = wLen/(radar.instrument_parameters['prt']['data'].copy()[0]*4)
            if vLow < 1.:
                tmp = vLow * 10**3
                if tmp < 1. or tmp > 100.:
                    print('\tYou need to check the PRT for DOW.')
                    print('\tProceeding with a PRT of {:.2f} us.'.format(
                            radar.instrument_parameters['prt']['data'].copy()[0]))
                    tmp = vLow
                vLow = tmp
        
        elif 'NOXP' in radar_file:
            vLow = wLen/(radar.instrument_parameters['prt']['data'].copy()[0] * 10**-3 *4)
            if vLow < 1.:
                tmp = vLow * 10**3
                if tmp < 1. or tmp > 100.:
                    print('\tYou need to check the PRT for NOXP.')
                    print('\tProceeding with a PRT of {:.2f} us.'.format(
                            radar.instrument_parameters['prt']['data'].copy()[0]))
                    tmp = vLow
                vLow = tmp
        else:
            vLow = wLen/(radar.instrument_parameters['prt']['data'].copy()[0]*4)
            if vLow > nyq:
                vLow = vLow /3.
            if vLow > nyq:
                print('\tFolding Nyquist cannot be determined.')
                break
            if vLow < 1.:
                tmp = vLow * 10**3
                if tmp < 1. or tmp > 100.:
                    print('\tYou need to check the PRT for NOXP.')
                    print('\tProceeding with a PRT of {:.2f} us.'.format(
                            radar.instrument_parameters['prt']['data'].copy()[0]))
                    tmp = vLow
                vLow = tmp
        
        ratioLow = nyq/vLow
        ratioLow = int(np.round(ratioLow,0))
        ratioHigh = ratioLow + 1
        
        vHigh = nyq/ratioHigh
        foldNyq = vHigh + vLow
        
        if stag:
            print('\tFolding Nyquist for dual-PRF corrections is {}'.format(foldNyq))
            if 'STP' not in radar.fields.keys():
                radar.add_field('STP',radar.fields[velOut].copy())
            radar.fields['STP']['data'] = prtCorrect(radar,velField = velOut,fnyq=foldNyq)
        
            despeckle(radar,vel_field=velOut,ngates=despk)
        
        if dual:
            print('\tLow and high Nyquists were calculated to be:  {:.3f}  {:.3f}'.format(vHigh,vLow))
            if useNew_dPRF:
                prfCorrect(radar,velField = velOut,nyqL=vLow,nyqH=vHigh)
            else:
                prfCorrectOld(radar,velField = velOut,nyqL=vLow,nyqH=vHigh)        
        
        

    #KDP Recalculation
    if kdpCalcFast or kdpCalcSlow:
        if str(phin) not in radar.fields.keys():
            print('No PHIDP field found. Skipping KDP calculation.')
        else: 
            try:
                if 'PHIC' not in radar.fields.keys():
                    radar.add_field('PHIC',radar.fields[phin].copy())
                if 'KDPC' not in radar.fields.keys():
                    radar.add_field('KDPC',radar.fields[kdpn].copy())
                if kdpCalcSlow:
                    radar.fields['PHIC']['data'],radar.fields['KDPC']['data']=\
                        kdpSlow(radar.fields[phin]['data'].copy(),
                                radar.fields[rhon]['data'].copy(),
                                radar.range['data'].copy(),3000.)
                if kdpCalcFast:
                    radar.fields['PHIC']['data'],radar.fields['KDPC']['data']=\
                        kdpFast(radar.fields[phin]['data'].copy(),
                                radar.fields[rhon]['data'].copy(),
                                radar.range['data'].copy(),3000.)
            except:
                print('Something went wrong with the KDP calculation(s). Skipping...')
                errorList.append(radar_file)
                np.savetxt('error.txt',np.array(errorList),fmt='%s')
                continue
    
    if attenCorr:
        if 'AC' not in radar.fields.keys():
            radar.add_field('AC',radar.fields[refOut].copy())
        if 'ACV' not in radar.fields.keys():
            radar.add_field('ACV',radar.fields[refOut].copy())
        if 'DZV' not in radar.fields.keys():
            radar.add_field('DZV',radar.fields[refOut].copy())
        if 'DZAC' not in radar.fields.keys():
            radar.add_field('DZAC',radar.fields[refOut].copy())
        if 'DZACV' not in radar.fields.keys():
            radar.add_field('DZACV',radar.fields[refOut].copy())
        if 'ZDAC' not in radar.fields.keys():
            radar.add_field('ZDAC',radar.fields[zdrn].copy())
        
        ac,acv,dzac,dzv,dzavc,zdac = calculate_attenuation(radar, 0, debug=True,
                                        doc=15, fzl=5000.0,rhv_min=0.5, ncp_min=0.0, 
                                        a_coef=0.06, beta=0.8,refl_field=refOut, ncp_field=snrn, 
                                        rhv_field=rhon,phidp_field='PHIC',zdr_field=zdrn, 
                                        spec_at_field='AC',corr_refl_field='DZAC',corr_zdr=True)
        
        radar.fields['AC']['data'],radar.fields['ACV']['data'],\
            radar.fields['DZAC']['data'],radar.fields['DZV']['data'],\
            radar.fields['DZACV']['data'],radar.fields['ZDAC']['data']=\
            ac['data'],acv['data'],dzac['data'],dzv,dzavc['data'],zdac
            
    
    if sfcRmv:
        print('- Running surface masking routine...')
        
        azPerSwp = 360
        
        bins = np.arange(azPerSwp)

        
        for iAf,(iaS,iaE,iaP) in enumerate(zip(applyAggStrtT,applyAggEndT,volStrmPos)):
            if filesDT[i] >= iaS and filesDT[i] <= iaE:
                tempStrmPos = iaP
                tempAggFn = outFnAgg[iAf]
                break
        
        print('\tImporting aggregate sweep file: {}'.format(os.path.split(tempAggFn)[1]))

        radAgg = pyart.io.read_cfradial(tempAggFn)
        aggMeanDBZ = radAgg.fields[refOut + '_mean']['data'].copy()
        aggStdDBZ = radAgg.fields[refOut + '_std']['data'].copy()
        
        
        if svExtra:
            radar.add_field_like(refOut,refOut + '_1A',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_1A',radar.fields[velOut]['data'].copy(),replace_existing=True)

        radar = airMaskSfc(radar,refOut,velOut,bins,aggMeanDBZ,aggStdDBZ,tempStrmPos,
                           aggMagThrsh=aggMagThrsh,sfcVel=sfcVel,
                           mskAzBuf_strm=mskAzBuf_strm,mskAzBuf_clr=mskAzBuf_clr,
                           mskDifThrsh=mskDifThrsh,altAdj=altAdj,zCsumThrsh=zCsumThrsh,
                           svExtra=svExtra)
                            
        if svExtra:
            radar.add_field_like(refOut,refOut + '_1B',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_1B',radar.fields[velOut]['data'].copy(),replace_existing=True)
                            
    if g2dFilt and not preDealThrsh:
        print('- Running radar noise filtering...')
        
        if svExtra:
            radar.add_field_like(refOut,refOut + '_2A',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_2A',radar.fields[velOut]['data'].copy(),replace_existing=True)
        
        radar = radFilt(radar,thrshVarID,refOut,velOut,thrshFillVal,
                        filtStdX=filtStdX,filtStdY=filtStdY,
                        smthThrsh=smthThrsh,retainShrtPulse=retainShrtPulse,
                        svExtra=svExtra)
        
        if svExtra:
            radar.add_field_like(refOut,refOut + '_2B',radar.fields[refOut]['data'].copy(),replace_existing=True)
            radar.add_field_like(velOut,velOut + '_2B',radar.fields[velOut]['data'].copy(),replace_existing=True)
    
    
    if clrVelChk:
        print('- Running clear air velocity checks/masking...')
        
        if svExtra:
            radar.add_field_like(velOut,velOut + '_3A',radar.fields[velOut]['data'].copy(),replace_existing=True)
        
        radar = caVelChk(radar,refOut,velOut,simVelID,refThrsh=refThrsh,
                         velDifThrsh=velDifThrsh,svExtra=svExtra)
        
        if svExtra:
            radar.add_field_like(velOut,velOut + '_3B',radar.fields[velOut]['data'].copy(),replace_existing=True)
    
    
    if rmvSpokes or any('pre' in rF for rF in sorted(glob(radar_folder+'/cfrad.*'))):
        radar_file_new = radar_file[:-11]+'_corr.nc'
    else:
        radar_file_new = radar_file[:-3]+'_corr.nc'
    
    date = radar_file_new.split('.')[1][:8]
    #pyart.io.write_cfradial(radar_file_new,radar,format='NETCDF4_CLASSIC',arm_time_variables=True)
    pyart.io.write_cfradial(radar_file_new,radar,arm_time_variables=True)
    print('\nCfRadial written...')

    if not deal4d_vol:
        if frmat == 'None':
            pass
        else:
            call(['RadxConvert -{} -outdir "{}" -f "{}" &>/dev/null'.format(frmat,radar_file[:radar_file.rindex('cf')-1],radar_file_new)],shell=True)
            #run('/usr/local/bin/RadxConvert -f {} -{} -outdir {} >/dev/null'.format(frmat,radar_folder,radar_file_new),shell=True)
            print ('{} written...'.format(frmat))
    if deal4d_vol: 
        call(['RadxConvert -{} -outdir "{}" -disag -f "{}"'.format(frmat,radar_file[:radar_file.rindex('cf')-1],radar_file_new)],shell=True)
        filePrev = radar_file_new
        
    
    timeF = dt.now()
    tDiff = timeF-time0
    
    print('DONE --> Seconds to edit sweep: {:.1f} . . . Estimated time remaining: {} '.format(tDiff.total_seconds(),str(tDiff*(len(radFiles)-(i+1)))[:-7]))    
        

wholeTimeF = dt.now()
wholeTDiff = wholeTimeF-wholeTime0

print('\n\nCIMMS-RadarQC complete. Total run time: {}'.format(wholeTDiff))