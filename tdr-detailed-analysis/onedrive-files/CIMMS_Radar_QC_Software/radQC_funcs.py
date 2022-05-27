########################
##   radQC_funcs      ##
##   VERSION: 1.3.8   ##
########################

from __future__ import division
import numpy as np
from warnings import warn
import pyart
from pyart.config import *
import pyart.correct.phase_proc as phase_proc
import matplotlib.pyplot as plt
import sys
import datetime
from astropy.convolution import Gaussian2DKernel, convolve
from scipy.optimize import curve_fit
from scipy.ndimage.filters import median_filter as mdfilt
from scipy.ndimage.filters import maximum_filter as mxfilter
from scipy.ndimage.filters import minimum_filter as mnfilter
from scipy.ndimage.filters import uniform_filter
from scipy.interpolate import UnivariateSpline
from scipy.stats import norm
from scipy.signal import savgol_filter
from scipy.integrate import cumtrapz
from scipy.ndimage.measurements import label


print("""

        ===================================================
        ==============     CIMMS-RadarQC     ==============
        ==============     Version 1.3.8     ==============
        ===================================================


""")


def despeckle(radar,vel_field=None, corr_vel_field=None, bad=-32768, ngates=4,AZ=False):
    """
    data = 1-D or 2-D array of radar data (e.g., reflectivity) to despeckle.
           If 1-D, should be a single ray. If 2-D, first dimension should be
           azimuth/elevation and second should be range.
    bad = Bad data value to check against
    ngates = Number of contiguous good data gates required along ray for no
             masking to occur
    Returns mask w/ same shape as data. True = bad.
    """ 
    if ngates == 0:
        return
    #vel_field, corr_vel_field = _parse_fields(vel_field, corr_vel_field)
    data = radar.fields[vel_field]['data'].filled(fill_value = -32768)
    #data = vdata.copy()     # dealiased velocities
    array = 1.0 * data
    mask = array == bad
    array[mask] = 0
    array[~mask] = 1
    if np.ndim(array) == 2:
        for az in np.arange(np.shape(array)[0]):
            mask[az, :] = _despeckle_ray(array[az, :], bad, ngates)
        if AZ==True:
            for r in np.arange(np.shape(array)[1]):
                mask[:, r] = _despeckle_ray(array[:, r], bad, int(ngates))
    elif np.ndim(array) == 1:
        mask = _despeckle_ray(array, bad, ngates)
    else:
        warn('Need 1- or 2-D array, #Failing ...')
        return
    data_mask = np.ma.masked_array(data,mask=(mask == True))
    radar.fields[vel_field]['data'] = np.ma.masked_array(data_mask,mask=(data_mask == -32768))
    #corr_vel = get_metadata(corr_vel_field)
    #corr_vel['data'] = data
    #return corr_vel
    return 


def _despeckle_ray(ray, bad, ngates):
    """
    Arguments similar to despeckle(). Should only be sent rays.
    """
    lab, nf = label(ray)
    if nf > 0:
        labels = lab[lab > 0]
        for feat in np.linspace(np.min(labels), np.max(labels),
                                num=np.max(labels)-np.min(labels)+1):
            featsz = np.size(labels[labels == feat])
            if featsz <= ngates:
                lab[lab == feat] = 0
    return lab == 0

#Author: Addison Alford
#For dual-PRF correction
def prfCorrectOld(radar,velField,nyqL=0,nyqH=0):
    timeo = datetime.datetime.now()

    #Define search width
    prfDiff = abs(nyqH-nyqL)/2.
    
    totalFixed,totalCaught = 0,0
    for sweep_slice in radar.iter_slice():
        #Define window length
        for window in [3,5,3]:
            #Define window orientation
            for center,label in zip(['opp',False, True],['Left','Right','Center']):
                #Define n = 1,2,3 or 2n = 2,4,6
                for d in [2.,4.,6.]:
                    fill=np.nan
                    velMask = radar.fields[velField]['data'][sweep_slice].mask
                    #Construct difference field
                    diff = radar.fields[velField]['data'][sweep_slice].filled(fill_value=fill) \
                        - median_roll(radar.fields[velField]['data'][sweep_slice].filled(fill_value=\
                        fill),window=window,center=center)
                    
                    #Fill masked points with 0 and save original difference field
                    diff = np.ma.masked_invalid(diff).filled(fill_value=0)
                    diffo = diff
                    
                    #Construct gaussian
                    mu,std = norm.fit(diffo)
                    thresh = 3*std
                    if thresh > nyqL:
                        thresh = 2*std
                    
                    #Mask inside mu +/- 3sigma
                    diffMask = np.ma.masked_inside(diffo,mu-thresh,mu+thresh)
                    
                    #Begin search for errors
                    div = np.ma.masked_where(diffMask.mask==True,diffo)/d
                    for nyq in [nyqL,nyqH]:
                        if nyq == min([nyqL,nyqH]):
                            divMask = np.ma.masked_outside(abs(div),(nyq)-prfDiff,(nyq)+prfDiff)
                        if nyq == max([nyqL,nyqH]):
                            divMask = np.ma.masked_outside(abs(div),(nyq)-prfDiff,(nyq)+prfDiff)
                        diff=np.ma.masked_where(divMask.mask==True,div)
                        
                        skipRad = 0
                        if diff.compressed().shape[0] ==0:
                            #print(('Skipping radial step for window = {}, '
                            #       'orientation = {}, and 2n = {:d}.').format(window,label,int(d)))
                            skipRad = 1
                        
                        if skipRad == 0:
                            #Try correcting errors of positive or negative magnitudes
                            #relative to surrounding points. If max-min values in
                            #window are not reduced, return to original value.                            
                            
                            #Sign of difference
                            sign = np.sign(diff)
                            pos = np.ma.masked_where(sign == -1.0,diff)
                            neg = np.ma.masked_where(sign == 1.0, diff)
                            vel = radar.fields[velField]['data'][sweep_slice].copy()
                            
                            #Check positive errors
                            mask = pos.mask == False

                            #Max-min before correction
                            stdBefore = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            
                            #Correction
                            vel[mask] = vel[mask] - d*nyq
                            
                            #Max-min after correction
                            stdAfter = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            
                            stdDiff = np.zeros_like(vel)
                            stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                            
                            #Where did correction increase the local max-min?
                            stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                            
                            signStd = np.sign(stdMask)
                            negStd = np.ma.masked_where(signStd == -1.0,signStd)
                            msk = negStd.mask == False
                            vel[msk] = vel[msk] + d*nyq
                            
                            #Compute statistics...
                            totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                            totalCaught = totalCaught + vel[msk].shape[0]
                            
                            #Repeat above for negative errors
                            mask = neg.mask == False
                            
                            stdBefore = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            
                            vel[mask] = vel[mask] + d*nyq
                            
                            stdAfter = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            
                            stdDiff=np.zeros_like(vel)
                            stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                            
                            stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                            
                            signStd = np.sign(stdMask)
                            negStd = np.ma.masked_where(signStd == -1.0,signStd)
                            msk = negStd.mask == False
                            vel[msk] = vel[msk] - d*nyq
                            
                            totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                            totalCaught = totalCaught + vel[msk].shape[0]
                            
                            #Update velocity field.
                            radar.fields[velField]['data'][sweep_slice] = np.ma.masked_where(velMask == True,vel)
                    
                    #Same as above, but for the azimuthal step.        
                    fill = np.nan
                    velMask = radar.fields[velField]['data'][sweep_slice].mask
                    
                    #Construct difference field
                    diff = radar.fields[velField]['data'][sweep_slice].filled(fill_value=fill) \
                        - median_rollaz(radar.fields[velField]['data'][sweep_slice].filled(fill_value=\
                        fill),window=window,center=center)
                    
                    diff = np.ma.masked_invalid(diff).filled(fill_value=0)
                    diffo = diff
                    
                    #New distribution
                    cdf = norm.cdf(diffo)
                    mu,std = norm.fit(diffo)
                    thresh = 3*std
                    if thresh > nyqL:
                        thresh = 2*std
                    
                    #Mask inside mu +/- 3sigma
                    diffMask = np.ma.masked_inside(diffo,mu-thresh,mu+thresh)
                    
                    #Proceed with correction at low and high Nyquists
                    div = np.ma.masked_where(diffMask.mask==True,diffo)/d
                    for nyq in [nyqL,nyqH]:
                        if nyq == min([nyqL,nyqH]):
                            divMask = np.ma.masked_outside(abs(div),(nyq)-prfDiff,(nyq)+prfDiff)
                        if nyq == max([nyqL,nyqH]):
                            divMask = np.ma.masked_outside(abs(div),(nyq)-prfDiff,(nyq)+prfDiff)
                        diff=np.ma.masked_where(divMask.mask==True,div)
                        
                        skipRad = 0
                        if diff.compressed().shape[0] ==0:
                            #print(('Skipping radial step for window = {}, '
                            #       'orientation = {}, and 2n = {:d}.').format(window,label,int(d)))
                            skipRad = 1
                        
                        if skipRad == 0:
                            sign = np.sign(diff)
                            pos = np.ma.masked_where(sign == -1.0,diff)
                            neg = np.ma.masked_where(sign == 1.0, diff)
                            vel = radar.fields[velField]['data'][sweep_slice].copy()
                            
                            #Negative points
                            mask = pos.mask == False
                            stdBefore = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            vel[mask] = vel[mask] - d*nyq
                            stdAfter = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            stdDiff = np.zeros_like(vel)
                            stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                            stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                            
                            signStd = np.sign(stdMask)
                            negStd = np.ma.masked_where(signStd == -1.0,signStd)
                            msk = negStd.mask == False
                            vel[msk] = vel[msk] + d*nyq
                            totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                            totalCaught = totalCaught + vel[msk].shape[0]
                            
                            #Positive points
                            mask = neg.mask == False
                            stdBefore = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            vel[mask] = vel[mask] + d*nyq
                            stdAfter = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                            stdDiff = np.zeros_like(vel)
                            stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                            stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                            
                            signStd = np.sign(stdMask)
                            negStd = np.ma.masked_where(signStd == -1.0,signStd)
                            msk = negStd.mask == False
                            vel[msk] = vel[msk] - d*nyq
                            
                            #Update stats
                            totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                            totalCaught = totalCaught + vel[msk].shape[0]
                            
                            #Update field
                            radar.fields[velField]['data'][sweep_slice] = np.ma.masked_where(velMask == True,vel)
    print('\tTotal Points Corrected = {:d}'.format(totalFixed))
    print('\tTotal Points Detected After False Correction = {:d}'.format(totalCaught))
    timef = datetime.datetime.now()

    print('\tSeconds to run dual-PRF algorithm {:.1f}'.format((timef-timeo).total_seconds()))
    return

def prfCorrect(radar,velField,nyqL=0,nyqH=0,maxNumIts = 1):
    timeo = datetime.datetime.now()
    pointsCaught = 0
    prfDiff = abs(nyqH-nyqL)/2
    for sweep_slice in radar.iter_slice():
        
        #Radial step first

        #Retain the original mask for the velocity field.
        maskO = radar.fields[velField]['data'][sweep_slice].mask.copy() #Added .copy. Fixed 23 Oct '19.
        
        #Test all intervals for each nyquist.

        #uniCount added as a way to define the maximum number of iterations.
        #maxNumIts is defined as an input now with a default of 1.
        uniCount=0
        while uniCount < maxNumIts:
            pointsCaughtPrev = pointsCaught
            for radFilter,azFilter,filterType in zip([51,21,21],[21,21,21],['both','radial','az']):
                
                #30 October 2019: The radial and azimuthal steps are now done before any correction is applied.
                #Both steps are performed indepently to retrieve the local error (smoothed radial or azimuthal
                #data are subtracted from radial or azimuthal linearly interpolated data to get local error).
                #Then the final error is the minimum error at each gate. The minimum assures both radial
                #and azimuthal identifications of outliers are indicating an erroeous point.
                #After the minimum step, the radial and azimuthal steps are then performed independtly to
                #capture radial or azimuthally oriented errors that still remain.

                ##### RADIAL STEP #####

                vel = radar.fields[velField]['data'][sweep_slice].copy()#Added .copy. Fixed 23 Oct '19.
                r = radar.range['data']
                
                #Create linearly interpolated sweeps in radial direction.
                velInterp = np.empty_like(vel.copy())
                for ray,num in zip(vel.copy(),range(velInterp.shape[0])):
                    if ray.compressed().shape[0] > 1:
                        rIn = np.ma.masked_where(ray.mask==True,r.copy())
                        rayOut = ray.copy() #Added .copy. Fixed 23 Oct '19.
                        rayOut = np.interp(r,rIn.compressed(),ray.compressed())
                        velInterp[num] = rayOut
                    else: velInterp[num] = ray.copy()

                #Allocate velocity and error arrays for each nyquist interval and d.
                velTest = np.empty((6,vel.shape[0],vel.shape[1]))
                error = np.ones((6,vel.shape[0],vel.shape[1]))*np.nan
                count = 0

                #Construct smooth profile on linearly-interpolated velocities.
                #Changed mode to interp because nearest was causing edge errors to 
                #be the "valid" data beyond the original data.
                velSmooth = savgol_filter(velInterp.copy(), radFilter, 3, mode='interp',axis=1)
                
                #Difference was incorrectly vel-velSmooth. Changed to velInterp - velSmooth. 27 October.
                #Difference also independent of nyq or d. Moved up. 27 October.
                #30 October: Save radial difference field.
                diffRad = np.ma.masked_where(maskO==True,velInterp.copy() - velSmooth.copy())
                diffRad = np.ma.masked_invalid(diffRad).filled(fill_value=0)
                mu,std = norm.fit(diffRad.flatten())
                threshRad = 1.*std

                ##### AZIMUTH STEP #####

                vel = radar.fields[velField]['data'][sweep_slice].T.copy()#Added .copy. Fixed 23 Oct '19.
                r = np.array(range(radar.azimuth['data'][sweep_slice].shape[0]))

                
                #Create linearly interpolated sweeps in azimuth direction.
                velInterp = np.empty_like(vel.copy())
                for ray,num in zip(vel.copy(),range(velInterp.shape[0])):
                    if ray.compressed().shape[0] > 1:
                        rIn = np.ma.masked_where(ray.mask==True,r.copy())
                        rayOut = ray.copy() #Added .copy. Fixed 23 Oct '19.
                        rayOut = np.interp(r,rIn.compressed(),ray.compressed())
                        velInterp[num] = rayOut
                    else: velInterp[num] = ray.copy()

                #Construct smooth profile on linearly-interpolated velocities.
                #Changed mode to interp because nearest was causing edge errors to 
                #be the "valid" data beyond the original data.
                velSmooth = savgol_filter(velInterp.copy(), azFilter, 3, mode='interp',axis=1)

                #Difference was incorrectly vel-velSmooth. Changed to velInterp - velSmooth. 27 October.
                #Difference also independent of nyq or d. Moved up. 27 October.
                #30 October: Save azimuthal difference field.
                diffAz = np.ma.masked_where(maskO.T==True,velInterp.copy() - velSmooth.copy())
                diffAz = np.ma.masked_invalid(diffAz).filled(fill_value=0)


                mu,std = norm.fit(diffRad)
                threshAz = 1.*std

                #Removed the threhsolding on the standard deviation for now.
                thresh = 0#np.nanmean([threshRad,threshAz])

                #Final difference field is min of radial and azimuthal difference fields or
                #independtly radial or azimuthal. 30 October.
                if filterType == 'both':diff = np.nanmin(np.array([diffRad,diffAz.T]),axis=0)
                if filterType == 'radial':diff = diffRad
                if filterType == 'az':diff = diffAz.T

                #Now make new copy of vel from original.
                vel = radar.fields[velField]['data'][sweep_slice].copy()

                #Loop over low and high nyquists and over d = 2, 4, and 6 for errors of
                # +/-2nV_N1 or +/-2nV_N2
                for nyq in [nyqL,nyqH]:
                    for d in [2,4,6]:
                        #velNew : holds original velocity information
                        #velToInterp : holds original velocity information with corrections masked out
                        velNew = vel.copy()
                        velToInterp = vel.copy()
                        
                        #Find positive errors.
                        positiveIndices = np.where((abs(diff) > thresh) & (diff != np.nan) & (diff >= 2*nyqL - 2*prfDiff) & (diff < 6*(nyqH + prfDiff)))

                        #Correct positive errors
                        velNew[positiveIndices] = velNew[positiveIndices] - d*nyq
                        pointsCaught = pointsCaught + (np.array(positiveIndices[0]).shape[0])/6
                        
                        #Remove errors and fill with NaNs.
                        velToInterp[positiveIndices] = np.nan

                        #print(velBefore[positiveIndices],velAfter[positiveIndices])
                        #sys.exit()
                        velToInterp.mask[positiveIndices] = True
                            
                        #Find negative errors
                        negativeIndices = np.where((abs(diff) > thresh) & (diff != np.nan) & (diff < -2*nyqL + 2*prfDiff) & (diff >= -6*(nyqH + prfDiff)))
                        
                        #Correct positive errors
                        velNew[negativeIndices] = velNew[negativeIndices] + d*nyq
                        pointsCaught = pointsCaught + (np.array(negativeIndices[0]).shape[0])/6
                        
                        #Remove errors and fill with NaNs.
                        velToInterp[negativeIndices] = np.nan
                        velToInterp.mask[negativeIndices] = True
                        
                        #New interopoation of sweep without errors present.
                        velInterpTestRad = np.empty_like(velToInterp.copy())
                        velInterpTestAz = np.empty_like(velToInterp.T.copy())

                        #Do interpolation of velocities with identified errors masked out.
                        r = radar.range['data']
                        for ray,num in zip(velToInterp.copy(),range(velInterp.shape[0])):
                            if ray.compressed().shape[0] > 1:
                                rIn = np.ma.masked_where(ray.mask==True,r.copy())
                                rayOut = ray.copy() #Added .copy. Fixed 23 Oct '19.
                                rayOut = np.interp(r,rIn.compressed(),ray.compressed())
                                velInterpTestRad[num] = rayOut
                            else: velInterpTestRad[num] = ray.copy()
                        
                        r = np.array(range(radar.azimuth['data'][sweep_slice].shape[0]))
                        for ray,num in zip(velToInterp.T.copy(),range(velInterp.shape[1])):
                            if ray.compressed().shape[0] > 1:
                                rIn = np.ma.masked_where(ray.mask==True,r.copy())
                                rayOut = ray.copy() #Added .copy. Fixed 23 Oct '19.
                                rayOut = np.interp(r,rIn.compressed(),ray.compressed())
                                velInterpTestAz[num] = rayOut
                            else: velInterpTestAz[num] = ray.copy()

                        #Create new smoothed field off linearly interpolated 
                        velSmoothTestRad = savgol_filter(velInterpTestRad.copy(), radFilter, 3, mode='interp',axis=1)
                        velSmoothTestAz = savgol_filter(velInterpTestAz.copy(), azFilter, 3, mode='interp',axis=1)

                        #Store corrected velocity field for 2nV_low or 2nV_high.
                        velTest[count] = velNew.copy()
                        
                        #Store error of field.
                        #Originally subtracted vel-velSmooth. Should be 
                        #velInterpTest-velSmooth. Changed 28 October.
                        errorRad = abs(velNew.filled(fill_value=np.nan).copy() - velSmoothTestRad.copy())
                        errorAz = abs(velNew.T.filled(fill_value=np.nan).copy() - velSmoothTestAz.copy())
                        
                        if filterType == 'both': error[count] = np.nanmin(np.array([errorRad,errorAz.T]),axis=0)
                        if filterType == 'radial':error[count] = errorRad
                        if filterType == 'az':error[count] = errorAz.T
                        
                        count+=1
                
                #Changed to velInterp-velSmooth on 28 October.
                #velTest[-1],error[-1] = vel,diff


                #Fill NaNs with infinity so their error will be ignored.
                error = np.ma.masked_invalid(error).filled(fill_value=np.inf)
                
                #Create "final" velocity field to be passed back to the radar object.
                velFinal = vel.copy()

                #Find minimum error along axis = 0.
                #Identify which correction applied produces the lowest local error
                #Any points that were never identified will retain their original value.
                
                #Indices where minimum error exists.
                minargs = np.nanargmin(error,axis=0)

                #Create index arrays for r, theta axes (axes = 1 and 2)
                rgrid,azgrid = np.meshgrid(range(vel.shape[0]),range(vel.shape[1]))
                
                #Make tuples out of indices
                #30 October: transposed arrays so rgrid,azgrid are of proper shape.
                #Had the effect of placing the wrong minargs index into the wrong
                #r,theta space.
                indices = tuple(np.array([minargs.flatten(),rgrid.T.flatten(),azgrid.T.flatten()]))

                #Boom. Final velocity field obtained
                #Switched order of reshape. This was a massive issue as the array was
                #being constructed improperly and giving incorrectly ordered data back
                #to the radar object. Then, in the subsequent azimuth step, everything
                #would be corrected, but with new errors introduced.
                #Bug fixed 23 October.
                velFinal = velTest[indices].reshape((vel.shape[0],vel.shape[1]))
                        
                #Update the radar field.
                radar.fields[velField]['data'][sweep_slice] = np.ma.masked_where(maskO==True,\
                    velFinal.copy())

                      
            if int(pointsCaught-pointsCaughtPrev) < 5 or uniCount == 0:
                print(pointsCaught,pointsCaughtPrev,uniCount)    
                print("\tTOTAL POINTS CAUGHT %d"%(int(pointsCaught)))
                radar.fields[velField]['data'][sweep_slice] = np.ma.masked_where(maskO==True,\
                    radar.fields[velField]['data'][sweep_slice]).copy()
                timef = datetime.datetime.now()

                print("\tSeconds to run dual-PRF algorithm %.1f" %(timef-timeo).total_seconds())

                return
            else:
                uniCount+=1
                print(pointsCaught,pointsCaughtPrev)
                print("\tPOINTS CAUGHT IN THIS ROUND %d"%(int(pointsCaught-pointsCaughtPrev)))
                
                

#For staggered-PRT
def prtCorrect(radar,velField,fnyq=0):
    totalFixed,totalCaught = 0,0
    #Define window length
    for window in [3,5,3]:
        #Define window orientation and search width (rng)
        for center,label,rng in zip(['opp',False, True],['Left','Right','Center'],[0.5,0.5,0.5]):
            fill = np.nan
            velMask = radar.fields[velField]['data'].mask
            for d in [1,2,3]:
                #Construct difference field for radial step.
                velMask = radar.fields[velField]['data'].mask
                diff = radar.fields[velField]['data'].filled(fill_value=fill) \
                    - median_roll(radar.fields[velField]['data'].filled(fill_value=\
                    fill),window=window,center=center)
                diff = np.ma.masked_invalid(diff)
                diffo = diff.filled(fill_value=0)
                
                #Construct normal distribution and retain points.
                mu,std = norm.fit(diffo)
                thresh = 3*std
                diffMask = np.ma.masked_inside(diffo,mu-thresh,mu+thresh)
                diffMask = np.ma.masked_outside(abs(diffMask),(d*fnyq) - (fnyq*rng),(d*fnyq) + (fnyq*rng))
                diff = np.ma.masked_where(diffMask.mask == True,diffo)
                
                #Check to see if there are any points to be corrected.
                skipRad = 0
                if diff.compressed().shape[0] == 0:
                    print(('Skipping radial step for window = {}, '
                           'orientation = {}, and n = {:d}.').format(window,label,int(d)))
                    skipRad = 1
                
                if skipRad == 0:
                    #Proceed with corrections
                    sign = np.sign(diff)
                    
                    #Error with positive difference
                    pos = np.ma.masked_where(sign == -1.0,diff)
                    #Error with negative difference
                    neg = np.ma.masked_where(sign == 1.0, diff)
                    vel = radar.fields[velField]['data'].copy()
                    
                    #Correct negative errors
                    mask = pos.mask == False

                    #Max-min before correction
                    stdBefore = max_roll(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    
                    #Correction
                    vel[mask] = vel[mask] - d*fnyq
                    
                    #Max-min after correction
                    stdAfter = max_roll(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    
                    #Check for points where max-min increased
                    stdDiff = np.zeros_like(vel)
                    stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                    stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                    
                    signStd = np.sign(stdMask)
                    negStd = np.ma.masked_where(signStd == -1.0,signStd)
                    msk = negStd.mask == False
                    
                    #Return any points that inceased local max-min to originals
                    vel[msk] = vel[msk] + d*fnyq
                    
                    #Update stats.
                    totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                    totalCaught = totalCaught + vel[msk].shape[0]
                    
                    #Proceed with positive values
                    mask = neg.mask == False
                    
                    #Max-min before correction
                    stdBefore = max_roll(vel,window=window,center=center) - \
                                            min_roll(vel,window=window,center=center)
                    
                    #Correction
                    vel[mask] = vel[mask] + d*fnyq
                    
                    #Max-min after correction
                    stdAfter = max_roll(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    
                    #Check for points where max-min increased
                    stdDiff = np.zeros_like(vel)
                    stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                    stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                    
                    signStd = np.sign(stdMask)
                    negStd = np.ma.masked_where(signStd == -1.0,signStd)
                    msk = negStd.mask == False
                    vel[msk] = vel[msk] - d*fnyq
                    
                    #Update stats
                    totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                    totalCaught = totalCaught + vel[msk].shape[0]
                    
                    #Update field
                    radar.fields[velField]['data'] = np.ma.masked_where(velMask == True,vel)
                
                #Construct difference field for azimuthal step
                diff = radar.fields[velField]['data'].filled(fill_value=fill) \
                    - median_rollaz(radar.fields[velField]['data'].filled(fill_value=\
                    fill),window=window,center=center)
                diff = np.ma.masked_invalid(diff)
                diffo = diff.filled(fill_value=0)
                
                #Construct distribution and mask inside mu +/- 3sigma
                cdf = norm.cdf(diffo)
                mu,std = norm.fit(diffo)
                thresh = 3*std
            
                diffMask = np.ma.masked_inside(diffo,-1*thresh,thresh)
                diffMask = np.ma.masked_outside(abs(diffMask),(d*fnyq) - (fnyq*rng),(d*fnyq) + (fnyq*rng))
                diff = np.ma.masked_where(diffMask.mask==True,diffo)
                
                #Check if there are points to be corrected.
                if diff.compressed().shape[0] ==0:
                    print(('Skipping radial step for window = {}, '
                           'orientation = {}, and n = {:d}.').format(window,label,int(d)))
                    skipRad = 1
                
                if skipRad == 0:
                    #Proceed with correction
                    sign = np.sign(diff)
                    pos = np.ma.masked_where(sign == -1.0,diff)
                    neg = np.ma.masked_where(sign == 1.0, diff)
                    vel = radar.fields[velField]['data'].copy()
                    
                    #Negative points
                    mask = pos.mask == False
                    
                    #Try correction as above
                    stdBefore = max_rollaz(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    
                    vel[mask] = vel[mask] - d*fnyq
                    
                    stdAfter = max_rollaz(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    
                    #Check correction is good   
                    stdDiff = np.zeros_like(vel)
                    stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                    stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                    
                    signStd = np.sign(stdMask)
                    negStd = np.ma.masked_where(signStd == -1.0,signStd)
                    msk = negStd.mask == False
                    vel[msk] = vel[msk] + d*fnyq
                    
                    #Update stats
                    totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                    totalCaught = totalCaught + vel[msk].shape[0]
                    
                    #Positive points
                    mask = neg.mask == False
                    stdBefore = max_rollaz(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    vel[mask] = vel[mask] + d*fnyq
                    stdAfter = max_rollaz(vel,window=window,center=center) - \
                                    min_roll(vel,window=window,center=center)
                    stdDiff = np.zeros_like(vel)
                    stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                    stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                    
                    signStd = np.sign(stdMask)
                    negStd = np.ma.masked_where(signStd == -1.0,signStd)
                    msk = negStd.mask == False
                    vel[msk] = vel[msk] - d*fnyq
                    
                    #Update stats
                    totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                    totalCaught = totalCaught + vel[msk].shape[0]
                    
                    #Update field
                    radar.fields[velField]['data'] = np.ma.masked_where(velMask == True,vel)
    
    print('Total Points Corrected = {:d}'.format(totalFixed))
    print('Total Points Detected After False Correction = {:d}'.format(totalCaught))
    
def prtCorrect2d(radar,velField,fnyq=0,thresh=0,its = 2):
    hists = []
    time = radar.time['units'].split('T')[1].split(':')
    its = 0
    totalFixed,totalCaught = 0,0
    count = 0
    for window in [3,5]:
        for center,label in zip(['opp'],['Left']):
            fill=np.nan
            velMask = radar.fields[velField]['data'].mask
            for d in [1,2,3]:
                velMask = radar.fields[velField]['data'].mask
                #diff = radar.fields[velField]['data'].filled(fill_value=fill) \
                 #   - median_roll(radar.fields[velField]['data'].filled(fill_value=\
                 #   fill),window=window,center=center)
                diff =  radar.fields[velField]['data'].filled(fill_value=fill) - \
                            mdfilt(radar.fields[velField]['data'].filled(fill_value=fill),
                                   size=(window,window),mode='constant',cval=np.nan)
                diff = np.ma.masked_invalid(diff)
                diffo = diff.filled(fill_value=0)
                mu,std = norm.fit(diffo)
                thresh = 3*std
                    
                #if count == 0 or count == 25:
                #    #########
                #    bins = range(-50,51,1)
                #    size=10
                #    #diffCheck = radar.fields[velField]['data'].filled(fill_value=fill) \
                #    #- median_roll(radar.fields[velField]['data'].filled(fill_value=\
                #    #fill),window=3,center=True)
                #    diffCheck=diff
                #    pdf = norm.fit(diffCheck.flatten())
                #    p = norm.pdf(bins,mu,std)
                #    fig = plt.figure(1)
                #    plt.clf()
                #    fig.set_size_inches(6,6)
                #    ax=fig.add_subplot(111)
                #    a=ax.hist(diffCheck.flatten(),bins,range=(np.nanmin(diffCheck.flatten()),np.nanmax(diffCheck.flatten())),alpha=0.5,histtype='bar',ec='black',align='mid',normed=True,log=True)
                #    b=ax.plot(bins,p,'k',linewidth=2)
                #    c = ax.plot([thresh,thresh],[1e-6,1e0],'r-',linewidth=2)
                #    c = ax.plot([-1*thresh,-1*thresh],[1e-6,1e0],'r-',linewidth=2)
                #    ax.set_xlim(min(bins),max(bins))
                #    ax.set_ylim(1e-6,1e0)
                #    plt.xticks(fontsize=size)
                #    plt.yticks(fontsize=size)
                #    ax.set_xlabel('Velocity Difference from Median (m s$\mathregular{^{-1}}$)',fontsize=size)
                #    #print label,window,d
                #    if count==0:
                #        ax.set_title('Before Correction | Centered, 3 Gate Radial Median'  ,fontsize=size)
                #        fig.savefig('%s%s%sUTC_before.png'%(time[0],time[1],time[2][0:2]),dpi=300)
                #        print 'No Corrections %.2f, %.2f...'%(mu,std)
                #    if count==25:
                #        ax.set_title('After Correction | Centered, 3 Gate Radial Median'  ,fontsize=size)
                #        fig.savefig('%s%s%sUTC_final.png'%(time[0],time[1],time[2][0:2]),dpi=300)
                #        print 'All Corrections %.2f, %.2f...'%(mu,std)
                #    plt.close()
                #    #########
                    
                
                
                count+=1
                diffMask = np.ma.masked_inside(diffo,-1*thresh,thresh)
                diffMask = np.ma.masked_outside(abs(diffMask),(d*fnyq) - (fnyq/2),(d*fnyq) + (fnyq/2))
                diff=np.ma.masked_where(diffMask.mask==True,diffo)
                skipRad = 0
                if diff.compressed().shape[0] ==0:
                    print('Skipping radial step for window = {}, iteration = {:d}.'.format(window,its))
                    skipRad = 1
                
                if skipRad == 0:
                    sign = np.sign(diff)
                    pos = np.ma.masked_where(sign == -1.0,diff)
                    neg = np.ma.masked_where(sign == 1.0, diff)
                    vel = radar.fields[velField]['data'].copy()
                    mask = pos.mask == False

                    ####
                    #stdBefore = max_roll(vel,window=window,center=center)-min_roll(vel,window=window,center=center)
                    stdBefore = mxfilt(vel,size=(window,window))-mnfilt(vel,size=(window,window))
                    ####
                    vel[mask] = vel[mask] - d*fnyq
                    
                    ####
                    #stdAfter = max_roll(vel,window=window,center=center)-min_roll(vel,window=window,center=center)
                    stdAfter = mxfilt(vel,size=(window,window))-mnfilt(vel,size=(window,window))
                    stdDiff = np.zeros_like(vel)
                    stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                    stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                    
                    signStd = np.sign(stdMask)
                    negStd = np.ma.masked_where(signStd == -1.0,signStd)
                    msk = negStd.mask == False
                    vel[msk] = vel[msk] + d*fnyq
                    totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                    totalCaught = totalCaught + vel[msk].shape[0]
                    
                    ####
                    
                    mask = neg.mask == False
                    ####
                    stdBefore = mxfilt(vel,size=(window,window))-mnfilt(vel,size=(window,window))
                    ####
                    vel[mask] = vel[mask] + d*fnyq
                    ####
                    ####
                    stdAfter = mxfilt(vel,size=(window,window))-mnfilt(vel,size=(window,window))
                    stdDiff = np.zeros_like(vel)
                    stdDiff[mask] = np.round(stdBefore-stdAfter,0)[mask]
                    stdMask = np.ma.masked_where(stdDiff>=0,stdDiff)
                    
                    signStd = np.sign(stdMask)
                    negStd = np.ma.masked_where(signStd == -1.0,signStd)
                    msk = negStd.mask == False
                    vel[msk] = vel[msk] - d*fnyq
                    totalFixed = totalFixed + (vel[mask].shape[0] - vel[msk].shape[0])
                    totalCaught = totalCaught + vel[msk].shape[0]
                    ####
                    radar.fields[velField]['data'] = np.ma.masked_where(velMask == True,vel)
    print('Total Points Corrected = {:d}'.format(totalFixed))
    print('Total Points Detected After False Correction = {:d}'.format(totalCaught))
    return

def std_filter(var,var_mask,ref,refThresh=10,mthresh=0.0,sthresh=6,window=5,snr=False):
    varn = var.copy()
    varnRad = var.copy()
    varnAz = var.copy()
    for i in range(0,len(var[:,0])):
        #varnRad[i] = pan.rolling_std(var[i],window,center=True)
        varnRad[i] = stdev_roll(var[i],window=window,center=True)
    
    for i in range(0,len(var[0,:])):
        #varnAz[:,i] = pan.rolling_std(var[:,i],window,center=True)
        varnAz[:,i] = stdev_roll(var[:,i],window=window,center=True)
    
    varn = np.nanmax(np.array([varnRad,varnAz]),axis=0)
    varn = np.ma.masked_array(varn,mask=np.isnan(varn)==True)
    varn = np.ma.masked_where(ref>refThresh,varn)
    
    if snr: 
        varn_mask = np.ma.masked_array(varn,mask=(var_mask!=0.0)).filled(fill_value=0)
    #varn_mask = np.ma.masked_array(varn,mask=(var_mask>mthresh)).filled(fill_value=0)
    #var_fil = np.ma.masked_array(var,mask=(varn_mask>sthresh))
    
    var_fil = np.ma.masked_array(var,mask=(varn.filled(fill_value=0)>sthresh))
    var_fil = np.ma.masked_array(var_fil,mask=(var_mask<mthresh))
    return var_fil

#Below are rolling median, mean, max, and min caluclations using scipy.

# (DS 10/18/19) - Added as some functions in here were still reliant upon
#   the pandas rolling_std() function which no longer exists in newer versions
# This was adapted from a blog post addressing rolling std. dev. and
#   seems relatively robust. The author modified a proposed method from StackOverflow.
#   Both links are included below:
#       https://nickc1.github.io/python,/matlab/2016/05/17/Standard-Deviation-(Filters)-in-Matlab-and-Python.html
#       https://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
def stdev_roll(var, window=5,center=True):
    varn = var.copy()
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    c1 = uniform_filter(varn, window, mode='nearest', origin=origin)
    c2 = uniform_filter(varn*varn, window, mode='nearest', origin=origin)
    return np.sqrt(c2 - c1*c1)

def median_roll(var,window=6,center=True):
    varn = var.copy()
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[:,0])):
        varn[i] = mdfilt(var[i],size=window,origin=origin,mode='nearest')

    return varn

def median_rollaz(var,window=3,center=True):
    varn = var.copy()
    
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[0])):
        varn[:,i] = mdfilt(var[:,i],size=window,origin=origin,mode='nearest')
    return varn

def max_rollaz(var,window=3,center=True):
    varn = var.copy()
    
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[0])):
        varn[:,i] = mxfilter(var[:,i],size=window,origin=origin,mode='nearest')
    return varn
    
def max_roll(var,window=3,center=True):
    varn = var.copy()
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[:,0])):
        varn[i] = mxfilter(var[i],size=window,origin=origin,mode='nearest')

    return varn
    
def max_rollaz(var,window=3,center=True):
    varn = var.copy()
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[0])):
        varn[:,i] = mxfilter(var[:,i],size=window,origin=origin,mode='nearest')

    return varn
    

def min_roll(var,window=3,center=True):
    varn = var.copy()
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[:,0])):
        varn[i] = mnfilter(var[i],size=window,origin=origin,mode='nearest')

    return varn
    
def min_rollaz(var,window=3,center=True):
    varn = var.copy()
    if center == True: 
        origin = 0 
    if center == False: 
        origin = int(window/2)
    if center == 'opp': 
        origin = int(window/2) *-1
    
    for i in range(0,len(var[0])):
        varn[:,i] = mnfilter(var[:,i],size=window,origin=origin,mode='nearest')

    return varn

def rad_interp(var):
    varn = var.copy()
    for i in range(0,len(var[:,0])):
        var_comp = var[i].compressed()
        if var_comp.shape[0]!=0:
            grid = np.ma.masked_array(range(0,len(var[i])),mask=var[i].mask==True).compressed()
            varn[i] = np.interp(range(0,len(var[i])),grid,var_comp)
    return varn

def grad_ndx(var,dea):
    var_grad = abs(np.gradient(var)[1]) + abs(np.gradient(var)[0])
    var_grad = np.ma.masked_array(var_grad,mask=np.isnan(var_grad)==True)
    var_grad = np.ma.masked_array(var_grad,dea!=0)
    var_grad = np.ma.masked_array(var_grad,var.mask==True)
    var_grad = np.ma.masked_array(var_grad,mask=var_grad==32768.0)
    vari = range(0,var.shape[0])
    varj = range(0,var.shape[1])
    varj,vari = np.meshgrid(varj,vari)
    vari = np.ma.masked_array(vari,var_grad.mask==True)
    vari = np.ma.masked_array(vari,var_grad<10)
    varj = np.ma.masked_array(varj,vari.mask==True)
    
    return np.array([vari.compressed(),varj.compressed()]).T


def var_rot_tel(var1,dea,deaF,nyq):
    varnew = var1.copy()
    var2 = var1-deaF*2*nyq
    ndx_list = grad_ndx(var1,dea)
    ndx_chng = []
    for window in [16,10,8,6,4]:
        indxS = window/2
        indxE = window/2 +1
        foot = np.ones((window-1))
        count = 1
        while count > 0:
            count = 0
            for ndxT in ndx_list:
                ndx = (ndxT[0],ndxT[1])
                if ndxT[1] > indxS and ndxT[1] < (var1.shape[1]-indxS) and ndx not in ndx_chng:
                    ndx = (ndxT[0],ndxT[1])
                    vartest1 = var1[ndx[0],ndx[1]-indxS:ndx[1]+indxE].copy()
                    vartest1[indxS] = var2[ndx].copy()
                    var1m = mdfilt(var1[ndx[0],ndx[1]-indxS:ndx[1]+indxE],footprint=foot)
                    var2m = mdfilt(vartest1,footprint=foot)
                    var1diff = abs(var1m-var1[ndx[0],ndx[1]-indxS:ndx[1]+indxE])[indxS]
                    var2diff = abs(var2m-vartest1)[indxS]
                    if var2diff < var1diff:
                        var1[ndx]=var2[ndx].copy()
                        ndx_chng.append(ndx)
                        count= count +1
            ndx_list = grad_ndx(var1,dea)
    return var1

def dealias_corr(var,deaF,nyq,wrap=True):
    foot1 = np.ones((25))
    foot2 = np.ones((8))
    foot3 = np.ones((3))
    varn = var.copy()+0.0
    for az in range(0,len(var[:,0])):
        ray = varn[az]+0.0
        dea = deaF[az]
        raymean = mdfilt(ray,footprint=foot1)
        raymean = np.ma.masked_array(raymean,mask=np.isnan(raymean)==True).filled(fill_value=0)
        raydiff = abs(ray-raymean)
        raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True).filled(fill_value=0)
        raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
        raydiff = np.ma.masked_array(raydiff,mask=np.max([abs(ray),abs(raymean)],axis=0)<(nyq*0.5))
        raym = np.ma.masked_array(ray,raydiff.mask==True)
        raym = np.ma.masked_array(raym,np.sign(raym)!=dea)
        rayg = np.ma.masked_array(ray,raym.mask==False).filled(fill_value=0)
        raynew = (raym-dea*nyq*2).filled(fill_value=0) + rayg
        #raynew = np.ma.masked_array(raynew,varn[az].mask==True)
        ray = raynew.copy()
        raydiff_old = raydiff.copy()
        raymean = mdfilt(ray,footprint=foot1)
        raymean = np.ma.masked_array(raymean,mask=np.isnan(raymean)==True).filled(fill_value=0)
        raydiff = abs(ray-raymean)
        raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True)
        raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
        raydiff = np.ma.masked_array(raydiff,mask=np.max([abs(ray),abs(raymean)],axis=0)<(nyq*0.5))
        while (raydiff.compressed().shape[0]<raydiff_old.compressed().shape[0]):
            raydiff_old = raydiff.copy()
            raym = np.ma.masked_array(ray,raydiff.mask==True)
            raym = np.ma.masked_array(raym,np.sign(raym)!=dea)
            rayg = np.ma.masked_array(ray,raym.mask==False).filled(fill_value=0)
            raynew = (raym-dea*nyq*2).filled(fill_value=0) + rayg
            #raynew = np.ma.masked_array(raynew,varn[az].mask==True)
            ray = raynew.copy()
            raymean = mdfilt(ray,footprint=foot1)
            raymean = np.ma.masked_array(raymean,mask=np.isnan(raymean)==True).filled(fill_value=0)
            raydiff = abs(ray-raymean)
            raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True)
            raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
            raydiff = np.ma.masked_array(raydiff,mask=np.max([abs(ray),abs(raymean)],axis=0)<(nyq*0.5))
    
        varn[az]=raynew
    azmed = mdfilt(varn,footprint=np.ones((3,1)))
    for az in range(0,len(var[:,0])):
        ray = varn[az]+0.0
        dea = deaF[az]
        raymeani = mdfilt(ray,footprint=foot2)
        raymeani = np.ma.masked_array(raymeani,mask=np.isnan(raymeani)==True).filled(fill_value=0)        
        raydiff = abs(azmed[az]-ray)
        raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True).filled(fill_value=0)
        raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
        raym = np.ma.masked_array(ray,raydiff.mask==True)
        rayg = np.ma.masked_array(ray,raym.mask==False).filled(fill_value=0)
        azsign = np.sign(azmed[az]-ray)
        rayt= (raym+azsign*nyq*2).filled(fill_value=0) + rayg
        #rayt= np.ma.masked_array(rayt,varn[az].mask==True)
        raymeanf = mdfilt(rayt,footprint=foot2)
        raymeanf = np.ma.masked_array(raymeanf,mask=np.isnan(raymeanf)==True).filled(fill_value=0)
        azsign = np.ma.masked_array(azsign,mask=(abs(ray-raymeani)<abs(rayt-raymeanf))).filled(fill_value=0)
        raynew = (raym+azsign*nyq*2).filled(fill_value=0) + rayg
        #raynew = np.ma.masked_array(raynew,varn[az].mask==True)
        ray = raynew.copy()
        varn[az]=raynew
    for az in range(0,len(var[:,0])):
        ray = varn[az]+0.0
        dea = deaF[az]
        raymean = mdfilt(ray,footprint=foot3)
        raymean = np.ma.masked_array(raymean,mask=np.isnan(raymean)==True).filled(fill_value=0)
        raydiff = abs(ray-raymean)
        raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True).filled(fill_value=0)
        raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
        raym = np.ma.masked_array(ray,raydiff.mask==True)
        raym = np.ma.masked_array(raym,np.sign(raym)!=dea)
        rayg = np.ma.masked_array(ray,raym.mask==False).filled(fill_value=0)
        raynew = (raym-dea*nyq*2).filled(fill_value=0) + rayg
        #raynew = np.ma.masked_array(raynew,varn[az].mask==True)
        ray = raynew.copy()
        raydiff_old = raydiff.copy()
        raymean = mdfilt(ray,footprint=foot3)
        raymean = np.ma.masked_array(raymean,mask=np.isnan(raymean)==True).filled(fill_value=0)
        raydiff = abs(ray-raymean)
        raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True)
        raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
        while (raydiff.compressed().shape[0]<raydiff_old.compressed().shape[0]):
            raydiff_old = raydiff.copy()
            raym = np.ma.masked_array(ray,raydiff.mask==True)
            raym = np.ma.masked_array(raym,np.sign(raym)!=dea)
            rayg = np.ma.masked_array(ray,raym.mask==False).filled(fill_value=0)
            raynew = (raym-dea*nyq*2).filled(fill_value=0) + rayg
            #raynew = np.ma.masked_array(raynew,varn[az].mask==True)
            ray = raynew.copy()
            raymean = mdfilt(ray,footprint=foot3)
            raymean = np.ma.masked_array(raymean,mask=np.isnan(raymean)==True).filled(fill_value=0)
            raydiff = abs(ray-raymean)
            raydiff = np.ma.masked_array(raydiff,mask=np.isnan(raydiff)==True)
            raydiff = np.ma.masked_array(raydiff,mask=raydiff<(nyq-1))
        varn[az]=raynew
    return varn

def line(x,m,b):
	return (m*x)+b    
	
def kdpSlow(var1,var2,rangeGate,filterLen):
    #Linear Least Squares fit using a 3 km window to fit PHI and KDP.
    #Threshold on RHOHV field
    dz = rangeGate[1]-rangeGate[0]
    medFiltLen = int(filterLen/dz)
    var = np.ma.masked_array(var1,mask=var2<0.3)
    print(var.shape)
    phiOut = np.ones_like(var)
    kdpOut = np.ones_like(var)
    deg=0
    for az in var:
        gate=0
        for r in rangeGate:
            if r > filterLen/2.:
                rMask=np.ma.masked_array(rangeGate,mask=az.mask==True)
                rMask = np.ma.masked_outside(rMask,r-filterLen/2.,r+filterLen/2.)
                tempVar = np.ma.masked_array(az,mask=rMask.mask==True).compressed()
                az = np.ma.masked_where(az==0.,az)
                if len(tempVar) > 3 and len(az[gate-4:gate].compressed()) > 3 and len(az[gate+1:gate+5].compressed()) > 3 and az.mask[gate] == False:
                    try:
                        coeff=curve_fit(line,rMask.compressed()/1000.,tempVar,p0=[1,1])[0]
                        phiOut[deg,gate] = line(r/1000.,coeff[0],coeff[1])
                        kdpOut[deg,gate] = coeff[0]/2.
                    except:
                        phiOut[deg,gate] = -32768.
                        kdpOut[deg,gate] = -32768.
                else:
                    phiOut[deg,gate] = -32768.
                    kdpOut[deg,gate] = -32768.
                    #print(coeff[0])
            else:
                phiOut[deg,gate] = -32768.
                kdpOut[deg,gate] = -32768.
            gate+=1
        deg+=1
    phiOut,kdpOut = np.ma.masked_array(phiOut,mask=phiOut==1.0),np.ma.masked_array(kdpOut,mask=kdpOut==1.0) 
    #phiOut = median_roll(phiOut,window=medFiltLen/4)
    return np.ma.masked_array(median_rollaz(phiOut,window=3),mask=phiOut==-32768.),np.ma.masked_array(median_rollaz(kdpOut),mask=kdpOut==-32768.)
    
def kdpFast(var1,var2,rangeGate,filterLen):
    #Spline interpolation of PHI field every filterLen m starting from gate = 0,1,2,3....
    #Takes median value of all spline interpolations to find a final solution
    #Runs a median filter on original phi (var1) field to remove Mie scattering
    #Spline is using a cubic smoother (k=3).
    #Threshold on RHOHV field
    var = np.ma.masked_array(var1,mask=var2<0.3)
    phiOut = np.ones_like(var)
    kdpOut = np.ones_like(var)
    dz = rangeGate[1]-rangeGate[0]
    medFiltLen = int(filterLen/dz)
    count=0
    
    for az in var:
        # (DS - 12/19/19) Switched to the local rolling median function from
        # the original code which relied on the now extinct Pandas rolling_median
        # **THIS HAS NOT YET BEEN TESTED**
        #azIn = pan.rolling_median(az,4)
        azIn = median_roll(az,window=4)
        azIn=np.ma.masked_array(azIn,mask=az.mask==True)
        #azIn=az
        splAZ=[]
        for i in range(0,medFiltLen):
            try:
                spl = UnivariateSpline(rangeGate[i:-1:medFiltLen],azIn[i:-1:medFiltLen],w=None,k=3)
                splAZ.append(spl(rangeGate))
            except: 
                return
            
        az_filt = np.median(np.array(splAZ),axis=0)
        phiOut[count]=np.ma.masked_array(az_filt,mask=az.mask==True)
        kdpOut[count,0:rangeGate.shape[0]-1]=np.diff(np.ma.masked_array(az_filt,mask=az.mask==True))/(2.*dz/1000.)
        count+=1
    return phiOut,kdpOut

            
            
        
def varsigndiff(var1,var2):
    deadiff = var1.filled(fill_value=0)-var2.filled(fill_value=0)
    dea = np.sign(deadiff)
    deas = np.ma.masked_array(dea,mask=dea==0.0)
    deaF = np.sign(rad_interp(deas))
    return dea, deaF
    
def var_fill(var):
    x_interp = np.ones((var.shape[1],))*range(0,var.shape[1])
    y_interp = np.ones((var.shape[0],))*range(0,var.shape[0])
    x_interp1,y_interp1 = np.meshgrid(x_interp,y_interp)
    x_interpm = np.ma.masked_array(x_interp1,mask=var.mask==True)
    var3 = var.copy()
    for i in range(0,len(var3[:,0])):
        try: 
            var3[i] = np.interp(x_interp,x_interpm[i].compressed(),var[i].compressed())
        except ValueError: 
            continue 
    varT = mdfilt(var3,size=4)
    return var.filled(fill_value = varT)
    
def calculate_attenuation(radar, z_offset, debug=False, doc=15, fzl=4000.0,
                          rhv_min=0.8, ncp_min=0.5, a_coef=0.06, beta=0.8,
                          refl_field=None, ncp_field=None, rhv_field=None,
                          phidp_field=None, zdr_field=None, spec_at_field=None,
                          corr_refl_field=None,corr_zdr=False):
    """
    Calculate the attenuation from a polarimetric radar using Z-PHI method.
    Parameters
    ----------
    radar : Radar
        Radar object to use for attenuation calculations.  Must have
        copol_coeff, norm_coherent_power, proc_dp_phase_shift,
        reflectivity_horizontal fields.
    z_offset : float
        Horizontal reflectivity offset in dBZ.
    debug : bool
        True to print debugging information, False supressed this printing.
    Returns
    -------
    spec_at : dict
        Field dictionary containing the specific attenuation.
    cor_z : dict
        Field dictionary containing the corrected reflectivity.
    Other Parameters
    ----------------
    doc : float
        Number of gates at the end of each ray to to remove from the
        calculation.
    fzl : float
        Freezing layer, gates above this point are not included in the
        correction.
    rhv_min : float
        Minimum copol_coeff value to consider valid.
    ncp_min : float
        Minimum norm_coherent_power to consider valid.
    a_coef : float
        A coefficient in attenuation calculation.
    beta : float
        Beta parameter in attenuation calculation.
    refl_field, ncp_field, rhv_field, phidp_field : str
        Field names within the radar object which represent the horizonal
        reflectivity, normal coherent power, the copolar coefficient, and the
        differential phase shift. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART
        configuration file.
    spec_at_field, corr_refl_field : str
        Names of the specific attenuation and the corrected
        reflectivity fields that will be used to fill in the metadata for
        the returned fields.  A value of None for any of these parameters
        will use the default field names as defined in the Py-ART
        configuration file.
    References
    ----------
    Gu et al. Polarimetric Attenuation Correction in Heavy Rain at C Band,
    JAMC, 2011, 50, 39-58.
    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phidp_field is None:
        # use corrrected_differential_phae or unfolded_differential_phase
        # fields if they are available, if not use differential_phase field
        phidp_field = get_field_name('corrected_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('unfolded_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('differential_phase')
    if spec_at_field is None:
        spec_at_field = get_field_name('specific_attenuation')
    if corr_refl_field is None:
        corr_refl_field = get_field_name('corrected_reflectivity')

    # extract fields and parameters from radar
    norm_coherent_power = radar.fields[ncp_field]['data']
    copol_coeff = radar.fields[rhv_field]['data']
    reflectivity_horizontal = radar.fields[refl_field]['data']
    proc_dp_phase_shift = radar.fields[phidp_field]['data']
    diff_ref = radar.fields[zdr_field]['data']
    nsweeps = int(radar.nsweeps)

    # determine where the reflectivity is valid, mask out bad locations.
    is_cor = copol_coeff > rhv_min
    is_coh = norm_coherent_power > ncp_min
    is_good = np.logical_and(is_cor, is_coh)
    mask = np.logical_not(is_good)
    refl = np.ma.masked_where(mask, reflectivity_horizontal + z_offset)

    # calculate initial reflectivity correction and gate spacing (in km)
    init_refl_correct = refl# + proc_dp_phase_shift * a_coef
    dr = (radar.range['data'][1] - radar.range['data'][0]) / 1000.0
    ZV = init_refl_correct - (diff_ref-1.1) 

    # create array to hold specific attenuation and attenuation
    specific_atten = np.zeros(reflectivity_horizontal.shape, dtype='float32')
    specific_atten_v = np.zeros(reflectivity_horizontal.shape, dtype='float32')
    atten = np.zeros(reflectivity_horizontal.shape, dtype='float32')
    atten_v = np.zeros(reflectivity_horizontal.shape, dtype='float32')
    atten_dp = np.zeros(reflectivity_horizontal.shape, dtype='float32')
    adp = np.zeros(reflectivity_horizontal.shape,dtype='float32')
    adp_final = np.zeros(reflectivity_horizontal.shape,dtype='float32')
    beta_estimate = np.zeros(reflectivity_horizontal.shape[0],dtype='float32')
    alphaList=[]
    dPhi = []

    for sweep in range(nsweeps):
        # loop over the sweeps
        if debug:
            print("Doing ", sweep)
        end_gate, start_ray, end_ray = phase_proc.det_process_range(radar, sweep, fzl, doc=doc)
        for i in range(start_ray, end_ray):
            # perform attenuation calculation on a single ray

            # extract the ray's phase shift and init. refl. correction
            ray_phase_shift = proc_dp_phase_shift[i, 0:end_gate]
            ray_init_refl = init_refl_correct[i, 0:end_gate]

            # perform calculation
            #last_six_good = np.where(is_good[i, 0:end_gate])[0][-6:]
            phidp_max = np.nanmax(ray_phase_shift)
            phidp_min = np.median(np.ma.masked_where(ray_phase_shift==-32768.,ray_phase_shift).compressed()[:10])
            delta_phi = phidp_max - phidp_min
            dPhi.append(delta_phi)
            #print phidp_max,phidp_min
            sm_refl = phase_proc.smooth_and_trim(ray_init_refl, window_len=5)
            reflectivity_linear = 10.0 ** (0.1 * beta * sm_refl)
            error=[]
            
            s_atten=np.zeros((end_gate,))
            a_range = np.arange(0.043,0.097,0.01)
            for a_coef in a_range:
                try:
                    self_cons_number = 10.0 ** (0.1 * beta * a_coef * delta_phi) - 1.0
                    I_indef = cumtrapz(0.46 * beta * dr * reflectivity_linear[::-1])
                    I_indef = np.append(I_indef, I_indef[-1])[::-1]
                    s_atten[0:end_gate-1]=cumtrapz(reflectivity_linear * self_cons_number /(I_indef[0] + (self_cons_number * I_indef)))*dr*2.0
                    phiMod = s_atten / a_coef
                    #print ray_phase_shift
                    error.append(np.nansum(abs(phiMod - ray_phase_shift)))
                except: 
                    error.append(np.nan)
            try:
                minIndex = error.index(np.nanmin(np.array(error)))
            except:
                minIndex = -999.
            #print a_range[minIndex]
            if minIndex != -999.:
                if delta_phi >= 30.:
                    alphaList.append(a_range[minIndex])
                    self_cons_number = 10.0 ** (0.1 * beta * a_range[minIndex] * delta_phi) - 1.0
                else:
                    alphaList.append(0.06)
                    self_cons_number = 10.0 ** (0.1 * beta * 0.06 * delta_phi) - 1.0
            else:
                self_cons_number = 0
                alphaList.append(np.nan)
            
            
            I_indef = cumtrapz(0.46 * beta * dr * reflectivity_linear[::-1])
            #print reflectivity_linear.shape
            #print ray_phase_shift.shape
            I_indef = np.append(I_indef, I_indef[-1])[::-1]
            # set the specific attenutation and attenuation
            specific_atten[i, 0:end_gate] = (reflectivity_linear * self_cons_number /(I_indef[0] + self_cons_number * I_indef))

            atten[i, :-1] = cumtrapz(specific_atten[i, :]) * dr * 2.0
            atten[i, -1] = atten[i, -2]
        #print alphaList
            
    # prepare output field dictionaries
    spec_at = get_metadata(spec_at_field)
    spec_at['data'] = atten
    spec_at['_FillValue'] = get_fillvalue()


    cor_z = get_metadata(corr_refl_field)
    cor_z['data'] = atten + reflectivity_horizontal + z_offset
    cor_z['data'].mask = init_refl_correct.mask
    cor_z['_FillValue'] = get_fillvalue()
    
    if corr_zdr == True:
        for sweep in range(nsweeps):
            # loop over the sweeps
            if debug:
                print("Doing ", sweep)
            end_gate, start_ray, end_ray = phase_proc.det_process_range(radar, sweep, fzl, doc=doc)
        
            for i in range(start_ray, end_ray):
                # perform attenuation calculation on a single ray

                # extract the ray's phase shift and init. refl. correction
                ray_phase_shift = proc_dp_phase_shift[i, 0:end_gate]
                ray_init_refl = ZV[i, 0:end_gate]

                # perform calculation
                #last_six_good = np.where(is_good[i, 0:end_gate])[0][-6:]
                phidp_max = np.max(ray_phase_shift)
                phidp_min = np.median(ray_phase_shift.compressed()[:10])
                delta_phi = phidp_max - phidp_min
                dPhi.append(delta_phi)
                #print phidp_max,phidp_min
                sm_refl = phase_proc.smooth_and_trim(ray_init_refl, window_len=5)
                reflectivity_linear = 10.0 ** (0.1 * beta * sm_refl)
                error=[]
                
                s_atten=np.zeros((end_gate,))
                a_range = np.arange(0.043,0.097,0.01)
                for a_coef in a_range:
                    try:
                        self_cons_number = 10.0 ** (0.1 * beta * a_coef * delta_phi) - 1.0
                        I_indef = cumtrapz(0.46 * beta * dr * reflectivity_linear[::-1])
                        I_indef = np.append(I_indef, I_indef[-1])[::-1]
                        s_atten[0:end_gate-1]=cumtrapz(reflectivity_linear * self_cons_number /(I_indef[0] + self_cons_number * I_indef))*dr*2.0
                        phiMod = s_atten / a_coef
                        #print ray_phase_shift
                        error.append(np.nansum(abs(phiMod - ray_phase_shift)))
                    except: 
                        error.append(np.nan)
                try:
                    minIndex = error.index(np.nanmin(np.array(error)))
                except:
                    minIndex = -999.
                #print a_range[minIndex]
                if minIndex != -999.:
                    if delta_phi >= 30.:
                        alphaList.append(a_range[minIndex])
                        self_cons_number = 10.0 ** (0.1 * beta * a_range[minIndex] * delta_phi) - 1.0
                    else:
                        alphaList.append(0.06)
                        self_cons_number = 10.0 ** (0.1 * beta * 0.06 * delta_phi) - 1.0
                else:
                    self_cons_number = 0
                    alphaList.append(np.nan)
                
                
                I_indef = cumtrapz(0.46 * beta * dr * reflectivity_linear[::-1])
                #print reflectivity_linear.shape
                #print ray_phase_shift.shape
                I_indef = np.append(I_indef, I_indef[-1])[::-1]
                # set the specific attenutation and attenuation
                specific_atten_v[i, 0:end_gate] = (reflectivity_linear * self_cons_number /(I_indef[0] + self_cons_number * I_indef))

                atten_v[i, :-1] = cumtrapz(specific_atten_v[i, :]) * dr * 2.0
                atten_v[i, -1] = atten_v[i, -2]
                
                adp = specific_atten[i,:]-specific_atten_v[i, :]
                
                atten_dp[i,:-1] = cumtrapz(np.ma.masked_where(adp<0,adp).filled(fill_value=0)) * dr * 2.0
                atten_dp[i,-1] = atten_dp[i,-2]
            #print alphaList
        
            
        # prepare output field dictionaries
        spec_at_v = get_metadata(spec_at_field)
        spec_at_v['data'] = atten_dp
        spec_at_v['_FillValue'] = get_fillvalue()


        cor_zv = get_metadata(corr_refl_field)
        cor_zv['data'] = atten_v + ZV + z_offset
        cor_zv['data'].mask = init_refl_correct.mask
        cor_zv['_FillValue'] = get_fillvalue()
        return spec_at, spec_at_v, cor_z, ZV, cor_zv, (diff_ref) + atten_dp 
    else:
        return spec_at, cor_z
        
def get_phidp_unf(radar, ncp_lev=0.4, rhohv_lev=0.6, debug=False, ncpts=20,\
          doc=-10, overide_sys_phase=False, sys_phase=-135,\
          nowrap=None, refl_field=None, ncp_field=None,\
          rhv_field=None, phidp_field=None):
    """
    Get Unfolded Phi differential phase. This function is pulled directly from
    pyart.correct.phase_proc with a number of modifications implemented.
    Parameters
    ----------
    radar : Radar
        The input radar.
    ncp_lev :
        Miminum normal coherent power level.  Regions below this value will
        not be included in the calculation.
    rhohv_lev :
        Miminum copolar coefficient level.  Regions below this value will not
        be included in the calculation.
    debug : bool, optioanl
        True to print debugging information, False to supress printing.
    ncpts : int
        Minimum number of points in a ray.  Regions within a ray smaller than
        this or beginning before this gate number are excluded from
        calculations.
    doc : int or None.
        Index of first gate not to include in field data, None include all.
    overide_sys_phase : bool, optional
        True to use `sys_phase` as the system phase. False will determine a
        value automatically.
    sys_phase : float, optional
        System phase, not used if overide_sys_phase is False.
    nowrap : or None
        Gate number where unwrapping should begin. `None` will unwrap all
        gates.
    refl_field ncp_field, rhv_field, phidp_field : str
        Field names within the radar object which represent the horizonal
        reflectivity, normal coherent power, the copolar coefficient, and the
        differential phase shift. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART
        configuration file.
    Returns
    -------
    cordata : array
        Unwrapped phi differential phase.
    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    if doc is not None:
        my_phidp = radar.fields[phidp_field]['data'][:, 0:doc]
        my_rhv = radar.fields[rhv_field]['data'][:, 0:doc]
        my_ncp = radar.fields[ncp_field]['data'][:, 0:doc]
        my_z = radar.fields[refl_field]['data'][:, 0:doc]
    else:
        my_phidp = radar.fields[phidp_field]['data']
        my_rhv = radar.fields[rhv_field]['data']
        my_ncp = radar.fields[ncp_field]['data']
        my_z = radar.fields[refl_field]['data']
    #t = time()
    if overide_sys_phase:
        system_zero = sys_phase
    else:
        system_zero = phase_proc.det_sys_phase(
            radar, ncp_field=ncp_field, rhv_field=rhv_field,
            phidp_field=phidp_field)
        if system_zero is None:
            system_zero = sys_phase
    cordata = np.zeros(my_rhv.shape, dtype=float)
    for radial in range(my_rhv.shape[0]):
        my_snr = phase_proc.snr(my_z[radial, :])
        notmeteo = np.logical_or(np.logical_or(
            my_ncp[radial, :] < ncp_lev,
            my_rhv[radial, :] < rhohv_lev), my_snr < 10.0)
        x_ma = ma.masked_where(notmeteo, my_phidp[radial, :])
        try:
            ma.notmasked_contiguous(x_ma)
            for slc in ma.notmasked_contiguous(x_ma):
                # so trying to get rid of clutter and small things that
                # should not add to phidp anyway
                if slc.stop - slc.start < ncpts or slc.start < ncpts:
                    x_ma.mask[slc.start - 1:slc.stop + 1] = True
            c = 0
        except TypeError:  # non sequence, no valid regions
            c = 1  # ie do nothing
            x_ma.mask[:] = True
        except AttributeError:
            # sys.stderr.write('No Valid Regions, ATTERR \n ')
            # sys.stderr.write(myfile.times['time_end'].isoformat() + '\n')
            # print x_ma
            # print x_ma.mask
            c = 1  # also do nothing
            x_ma.mask = True
        if 'nowrap' is not None:
            # Start the unfolding a bit later in order to avoid false
            # jumps based on clutter
            unwrapped = copy.deepcopy(x_ma)
            end_unwrap = phase_proc.unwrap_masked(x_ma[nowrap::], centered=False)
            unwrapped[nowrap::] = end_unwrap
        else:
            unwrapped = phase_proc.unwrap_masked(x_ma, centered=False)
        #end so no clutter expected
        system_max = unwrapped[np.where(np.logical_not(
            notmeteo))][-10:-1].mean() - system_zero
        unwrapped_fixed = np.zeros(len(x_ma), dtype=float)
        based = unwrapped-system_zero
        based[0] = 0.0
        notmeteo[0] = False
        based[-1] = system_max
        notmeteo[-1] = False
        '''#unwrapped_fixed[np.where(np.logical_not(based.mask))[0]] = \
            based[np.where(np.logical_not(based.mask))[0]]
        #if len(based[np.where(np.logical_not(based.mask))[0]]) > 11:
            unwrapped_fixed[np.where(based.mask)[0]] = \
                np.interp(np.where(based.mask)[0],
                          np.where(np.logical_not(based.mask))[0],
                          phase_proc.smooth_and_trim(based[np.where(
                              np.logical_not(based.mask))[0]]))
        #else:
            unwrapped_fixed[np.where(based.mask)[0]] = \
                np.interp(np.where(based.mask)[0],
                          np.where(np.logical_not(based.mask))[0],
                          based[np.where(np.logical_not(based.mask))[0]])
                          '''
        unwrapped_fixed=unwrapped
        if c != 1:
            cordata[radial, :] = unwrapped_fixed
        else:
            cordata[radial, :] = np.zeros(my_rhv.shape[1])
    if debug:
        print("Exec time: ", time() - t)
    return cordata


def rayZCsum(dbzArr):
    """
    This function calculates the cumulative sum of reflectivity in terms
    of Z (divided by 1x10^6 to ease interpretation) along each radial.    
    
    
    Parameters
    ----------
    dbzArr : float array
        Array containing the reflectivity values on which to perform
        the cumulative sum. If this array is coming from a Py-ART radar
        object and has not been modified otherwise, it is expected to be
        a numpy masked array.
    
    
    Returns
    -------
    zCsum : float array
        Array containing the cumulative sum of Z (divided by 1x10^6)
        along each ray.
    """

    DBZ = dbzArr.copy()
    
    # Undo the log
    z = 10. ** (DBZ / 10.)
    
    # Calculate the cumulative sum at each gate along each radial 
    # and then divide by 1x10^6 to return values easier to work with 
    zCsum = np.ma.cumsum(z,axis=1)/1000000.
    
    return zCsum




def rmvSpoke(radar,dbzID,velID,air,spokeDiff=10,maskRingOut=None,maskRingIn=None,
             dualPRTofst=False,svExtra=False):
    """
    This function is intended to identify and remove "spokes" from a radar sweep. Spokes
    are typically seen as groups of locally higher reflectivity along a radial (relative
    to the radials before and after). Gates are analyzed along each radial by considering
    the sign and magnitude of the difference between the target gate and the gates
    immediately before and after in azimuth. If the target gate is greater than both of
    these adjacent gates by at least the value defined by `spokeDiff`, or if the adjacent
    gates are both masked/invalid, then the target gate is flagged as a spoke member.
    Currently, this function does not consider the number of consecutive potential spoke
    members, and thus also tends to act as a pseudo-defreckling/despeckling method. Future
    versions will likely adjust this behavior to prevent removal of good echoes.
    
    
    Parameters
    ----------
    radar : Py-ART radar object
        Radar object containing the reflectivity and velocity fields to be processed.
    dbzID : string
        Name of the reflectivity variable to be processed.
    velID : string
        Name of the velocity variable to be processed.
    spokeDiff : float, optional
        Threshold difference (dBZ) between a target gate and the one preceding it 
        (in azimuth) at/beyond which the target gate is flagged as a spoke member.
        Default value is 10, which seems to work well with at least the TDR data 
        from TORUS19.
    maskRingOut : float, optional
        Range (in meters) from the radar at and beyond which all gates will be masked. 
        (e.g., a value of 45000 will mask all gates at ranges of 45000 meters and beyond)
        Useful for removing the often erroneous data near the max range of the radar. 
        A value of 0 or `None` (default) will turn off this feature.
    maskRingIn : float, optional
        Range (in meters) from the radar at and within which all gates will be masked.
        (e.g., a value of 5000 will mask all gates between the radar and 5000 meters range)
        Useful for removing the often erroneous data in gates nearest the radar. 
        A value of 0 or `None` (default) will turn off this feature.
    dualPRTofst : bool, optional
        Apply an offset correction to account for differences in reflectivity sensitivity 
        due to use of both a short and long pulse. Applies specifically to TORUS TDR data 
        where a +9.8 dBZ offset is added to reflectivity data within 450 m and 5625 m. 
        May be made more configurable by the user if/when other use cases are identified.
    svExtra : bool, optional
        Flag to enable output of intermediate/additional variables useful in assessing
        the performance of this editing step.
    
    
    Returns
    -------
    radar : Py-ART radar object
        Py-ART radar object containing original and QC'd data. This radar object can then be used for
        additional analysis or plotting, as well as output to finalized cfRadial files.
    """
    time0 = datetime.datetime.now()
    
    if air:
        az_raw = radar.rotation['data'].data
        roll = radar.roll['data'].data
        azmth = az_raw + roll
        azmth[azmth < 0] += 360
        azmth[azmth > 360] -= 360
    else:
        azmth = radar.azimuth['data'].data
        azmth[azmth < 0] += 360
        azmth[azmth > 360] -= 360
        
    rng = radar.range['data'].data

    DBZ = radar.fields[dbzID]['data'].copy()
    VEL = radar.fields[velID]['data'].copy()
    
    if dualPRTofst:
        DBZ[:,6:76] += 9.8
        
    if maskRingOut is not None and maskRingOut != 0:
        DBZ[:,rng >= maskRingOut] = np.ma.masked
        VEL[:,rng >= maskRingOut] = np.ma.masked
        
    if maskRingIn is not None and maskRingIn != 0:
        DBZ[:,rng <= maskRingIn] = np.ma.masked
        VEL[:,rng <= maskRingIn] = np.ma.masked
    
    spkMsk = np.full_like(DBZ.mask,False)
    prvPosDif = np.ones_like(DBZ.data)*np.nan

    DBZfld = DBZ.filled(fill_value=np.nan)

    #print('DBZ shape: {}'.format(DBZ.shape))
    #print('Mask shape: {}'.format(spkMsk.shape))

    timeIn0 = datetime.datetime.now()
    for ia in range(0,len(azmth)):
        if ia == 0:
            iaB = len(azmth) - 1
            iaE = ia + 1
        elif ia == len(azmth) - 1:
            iaB = ia - 1
            iaE = 0
        else:
            iaB = ia - 1
            iaE = ia + 1
        slc = [iaB,ia,iaE]
        for ir in range(0,len(rng)):
            tmpSlc = DBZfld[slc,ir]
        
            # If our slice is completely masked or the middle (target) value
            # is masked, we skip to the next slice. If the values to either
            # side of the target are masked, we go ahead and count the target
            # a potential spoke.
            if np.all(np.isnan(tmpSlc)) or np.isnan(tmpSlc[1]):
                continue
            elif np.all(np.isnan(tmpSlc[[0,2]])):
                spkMsk[ia,ir] = True
                continue
            #elif np.isnan(tmpSlc[0]) or np.isnan(tmpSlc[2]):
            #    continue
            else:
                if np.isnan(tmpSlc[0]):
                    tmpSlc[0] = tmpSlc[2]
                elif np.isnan(tmpSlc[2]):
                    tmpSlc[2] = tmpSlc[0]
                
                tmpDif = np.diff(tmpSlc)
                if np.array_equal(np.sign(tmpDif),[1,-1]):
                    prvPosDif[ia,ir] = tmpDif[0]
                    try: #
                        if tmpDif[0] >= spokeDiff:
                            spkMsk[ia,ir] = True
                    except:
                        print('Failed -- ia = {}  ir = {}'.format(ia,ir))
                        sys.exit()
    
    timeInF = datetime.datetime.now()
    tInDiff = timeInF-timeIn0
    
    #print('\t\tSeconds to loop: {:.1f} '.format(tInDiff.total_seconds())) 
    
    prvPosD_azCnt = np.nansum(prvPosDif,axis=1)
    
    DBZ_new = np.ma.masked_array(data=DBZ.data,mask=np.logical_or(spkMsk,DBZ.mask),copy=True)
    #VEL_new = np.ma.masked_array(data=VEL.data,mask=np.logical_or(spkMsk,VEL.mask),copy=True)
    VEL_new = np.ma.masked_array(data=VEL.data,mask=np.logical_or(DBZ_new.mask,VEL.mask),copy=True)
    radar.add_field_like(dbzID,dbzID,DBZ_new,replace_existing=True)
    radar.add_field_like(velID,velID,VEL_new,replace_existing=True)
    
    if svExtra:
        radar.add_field_like(dbzID,dbzID + '_0spkD',prvPosDif,replace_existing=True)
        dDBZ_spok = DBZ.filled(fill_value=0)-DBZ_new.filled(fill_value=0)
        radar.add_field_like(dbzID,dbzID + '_0sDd',dDBZ_spok,replace_existing=True)
    
    timeF = datetime.datetime.now()
    tDiff = timeF-time0
    
    #print('\t\tTotal seconds to remove spokes: {:.1f} '.format(tDiff.total_seconds())) 
    
    return radar
    
    


def genAgg(inFiles,strmPos,velOut,refOut,velIn,refIn):
    """
    This function combines a volume of radar data (currently only reflectivity, radial 
    velocity, and altitude are combined) into a single 3-dimensional array of shape 
    (numSweeps (`len(inFiles)`), numAzimuths (360), `gatesPerRay`). Not every sweep has 
    rays at exactly the same azimuths (rotation angles), so data from each sweep are   
    sorted into 1° bins spanning 0°-359°, and some number of gates within each 1° 
    radial (radar-dependent). Following aggregation of all volume sweeps into a common 
    array, the mean and standard deviation of reflectivity and radial velocity are 
    calculated at each gate within each 1° radial, yielding output arrays of 
    shape (numAzimuths (360), `gatesPerRay`). The primary purpose of this function is to 
    generate mean and standard deviation fields representative of radar observations in 
    clear air within a given volume. As such, the storm position (`strmPos`) string
    specifies whether the target storm was to the left ('L') or right ('R') of the
    aircraft heading. Data from the side opposite of `strmPos` is assumed to be
    representative of clear air conditions within the volume, with the mean and standard
    deviation data within the clear air side (0°-180° azimuth/rotation angle if `strmPos` 
    is 'L'; 180°-359° azimuth/rotation angle if `strmPos` is 'R') mirrored to the storm
    side. The result should be representative of the mean and standard deviation of
    reflectivity and radial velocity within the given volume had the radar been scanning
    in entirely clear air.
    The resulting aggregate stats fields are saved to a Py-ART radar object. Currently, 
    this step searches for a sweep nearest the midpoint of the volume with exactly 360
    azimuth angles, creates a radar object from that sweep and overwrites the data within
    with the aggregate sweep variables. This ensures the radar object has all the 
    necessary metadata used when generating an output CfRadial file. It is important to
    note that certain metadata of this output file will be incorrect (e.g., altitude, time),
    and will instead reflect the parameters of the original sweep. This is not an issue for
    most applications of the output, save for perhaps plotting functions. Future versions 
    of this function will produce a custom radar object where all metadata are accurate.
    
    
    Parameters
    ----------
    inFiles : list of strings
        List containing the paths to each input CfRadial file to be included in the
        calculation of the volume statistics.
    strmPos : single character string
        Side of the aircraft heading where the target storm/feature is located. Data from
        gates opposite this will be mirrored in the final output variables. Valid options
        include 'L' when the target storm/feature is to the left of the aircraft, and 
        'R' when the target storm/feature is to the right of the aircraft.
    velOut : string
        Name of the output velocity field. If it does not exist in the input file,
        the program will attempt to use `velIn`.
    refOut : string
        Name of the output reflectivity field. If it does not exist in the input file,
        the program will attempt to use `refIn`.
    velIn : string
        Name of the input velocity field. If it does not exist in the input file,
        the script will exit with an error.
    refIn : string
        Name of the input reflectivity field. If it does not exist in the input file,
        the script will exit with an error.
    
    Returns
    -------
    radar : Py-ART radar object
        Py-ART radar object with volume-aggregate statistical fields.
    """
    
    azPerSwp = 360
    
    bins = np.arange(azPerSwp)
    
    # Read in one of the input files to determine how many gates per ray there are and
    #   to see if the reflectivity and velocity variables exist
    tmpRadar0 = pyart.io.read_cfradial(inFiles[0])
    gatesPerRay = len(tmpRadar0.range['data'].data)
    
    swpAngl = tmpRadar0.fixed_angle['data'].data[0]
    
    
    if velOut in tmpRadar0.fields.keys():
        velID = velOut
    elif velIn in tmpRadar0.fields.keys():
        velID = velIn
    else:
        print('Neither {} nor {} exist (within at least the first CfRadial file). Exiting...'.format(velOut,velIn))
        sys.exit()
    
    if refOut in tmpRadar0.fields.keys():
        dbzID = refOut
    elif refIn in tmpRadar0.fields.keys():
        dbzID = refIn
    else:
        print('Neither {} nor {} exist (within at least the first CfRadial file). Exiting...'.format(refOut,refIn))
        sys.exit()
    
    del tmpRadar0
    
    aggDBZ = np.ones((len(inFiles),azPerSwp,gatesPerRay),dtype=np.float32)*np.nan
    aggVEL = np.ones((len(inFiles),azPerSwp,gatesPerRay),dtype=np.float32)*np.nan
    aggALT = np.ones((len(inFiles),),dtype=np.float32)*np.nan

    #fullSwpFIx = []

    aggMeanDBZ = np.full((azPerSwp,gatesPerRay),np.nan)
    aggMeanVEL = np.full((azPerSwp,gatesPerRay),np.nan)
    aggStdDBZ = np.full((azPerSwp,gatesPerRay),np.nan)
    aggStdVEL = np.full((azPerSwp,gatesPerRay),np.nan)

    if strmPos == 'L':
        posStr = 'left'
    elif strmPos == 'R':
        posStr = 'right'
    
    posRcnt = 0
    posLcnt = 0
    
    print('\tAggregating sweeps to generate statistical masking fields. Storm is to the {} of the radar.'.format(posStr))
    for i,radF in enumerate(inFiles):
        #if not i%10:
        print('\tSweep {} of {}...'.format(i+1,len(inFiles)))
        tmpRadar = pyart.io.read_cfradial(radF)
    
        tmpAz = tmpRadar.rotation['data'].data + tmpRadar.roll['data'].data
        tmpAz[tmpAz < 0] += 360
        tmpAz[tmpAz >= 360] -= 360

        #if len(tmpAz) == azPerSwp:
        #    fullSwpFIx.append(i)

        tmpDBZ = tmpRadar.fields[dbzID]['data'].copy().filled(fill_value=np.nan)
        tmpVEL = tmpRadar.fields[velID]['data'].copy().filled(fill_value=np.nan)

        tmpNewIx = (np.digitize(tmpAz,bins,right=False))-1

        aggDBZ[i,tmpNewIx,:] = tmpDBZ
        aggVEL[i,tmpNewIx,:] = tmpVEL
        aggALT[i] = tmpRadar.altitude_agl['data'][0]
        
        tmpMskCntL = np.ma.count_masked(np.ma.masked_invalid(aggDBZ[i,180:,:]))
        tmpMskCntR = np.ma.count_masked(np.ma.masked_invalid(aggDBZ[i,0:180:,:]))
        
        if tmpMskCntL > tmpMskCntR:
            posRcnt += 1
           #print('\t\tLeft side count greater than right by {} points'.format(tmpMskCntL-tmpMskCntR))
        if tmpMskCntR > tmpMskCntL:
            posLcnt += 1
           #print('\t\tRight side count greater than left by {} points'.format(tmpMskCntR-tmpMskCntL))
        

    if posRcnt > posLcnt:
        autoStrmPos = 'R'
    else:
        autoStrmPos = 'L'
        
    if strmPos != autoStrmPos:
        print('\t*!*!* User-specified storm position ( {} ) does not agree\n'
              '\t*!*!* with the automatically determined position ( {} )'.format(strmPos.upper(),autoStrmPos.upper()))

    aggDBZ = np.ma.masked_invalid(aggDBZ)
    aggVEL = np.ma.masked_invalid(aggVEL)
    aggALT = np.ma.masked_invalid(aggALT)

    aggMeanDBZ1 = np.ma.mean(aggDBZ,axis=0)
    aggMeanVEL1 = np.ma.mean(aggVEL,axis=0)
    aggStdDBZ1 = np.ma.std(aggDBZ,axis=0)
    aggStdVEL1 = np.ma.std(aggVEL,axis=0)

    if strmPos == 'L':
        aggMeanDBZ[0:180,:] = aggMeanDBZ1[0:180,:]
        aggMeanVEL[0:180,:] = aggMeanVEL1[0:180,:]
        aggStdDBZ[0:180,:] = aggStdDBZ1[0:180,:]
        aggStdVEL[0:180,:] = aggStdVEL1[0:180,:]
        if len(strmPos) == 1:
            aggMeanDBZ[180:,:] = np.flipud(aggMeanDBZ1[0:180,:])
            aggMeanVEL[180:,:] = np.flipud(aggMeanVEL1[0:180,:])
            aggStdDBZ[180:,:] = np.flipud(aggStdDBZ1[0:180,:])
            aggStdVEL[180:,:] = np.flipud(aggStdVEL1[0:180,:])

    if strmPos == 'R':
        aggMeanDBZ[180:,:] = aggMeanDBZ1[180:,:]
        aggMeanVEL[180:,:] = aggMeanVEL1[180:,:]
        aggStdDBZ[180:,:] = aggStdDBZ1[180:,:]
        aggStdVEL[180:,:] = aggStdVEL1[180:,:]
        if len(strmPos) == 1:
            aggMeanDBZ[0:180,:] = np.flipud(aggMeanDBZ1[180:,:])
            aggMeanVEL[0:180,:] = np.flipud(aggMeanVEL1[180:,:])
            aggStdDBZ[0:180,:] = np.flipud(aggStdDBZ1[180:,:])
            aggStdVEL[0:180,:] = np.flipud(aggStdVEL1[180:,:])


    radarAgg = initAggRad(gatesPerRay,azPerSwp,bins,swpAngl)

    radarAgg.add_field(dbzID + '_mean',initAggField(aggMeanDBZ,'dBZ'))
    radarAgg.add_field(velID + '_mean',initAggField(aggMeanVEL,'m/s'))
    radarAgg.add_field(dbzID + '_std',initAggField(aggStdDBZ,'dBZ'))
    radarAgg.add_field(velID + '_std',initAggField(aggStdVEL,'m/s'))

    #midPtFullSwp = fullSwpFIx[(np.abs(np.asarray(fullSwpFIx) - len(inFiles)/2)).argmin()]
    #radarAgg = pyart.io.read_cfradial(inFiles[midPtFullSwp],include_fields=[dbzID,velID]) # Chosen as it has same number of radials (360) as aggregate PDD fields
    #radarAgg.rotation['data'] = bins
    #radarAgg.roll['data'] = np.ma.masked_array(np.zeros((azPerSwp,)))

    #radarAgg.add_field_like(dbzID,dbzID + '_mean',aggMeanDBZ,replace_existing=True)
    #radarAgg.add_field_like(velID,velID + '_mean',aggMeanVEL,replace_existing=True)
    #radarAgg.add_field_like(dbzID,dbzID + '_std',aggStdDBZ,replace_existing=True)
    #radarAgg.add_field_like(velID,velID + '_std',aggStdVEL,replace_existing=True)

    ## Delete these fields as they were only needed to provide metadata for our aggregate fields
    #del radarAgg.fields[dbzID]
    #del radarAgg.fields[velID]
        
    
    return radarAgg


def initAggRad(numGates,numRays,azimuth,fixedAngle,numSweeps=1):
    aggRadar = pyart.testing.sample_objects.make_empty_ppi_radar(numGates,numRays,numSweeps)
    aggRadar.fixed_angle['data'][0] = fixedAngle
    aggRadar.azimuth['data'] = azimuth
    
    return aggRadar
    
def initAggField(dataArr,units,fillVal=-32768.0):
    outDict = {'units': units,
               '_FillValue': fillVal,
               'sampling_ratio': 1.0,
               'grid_mapping': 'grid_mapping',
               'coordinates': 'time range',
               'data': dataArr}
    
    return outDict


def airMaskSfc(radar,dbzID,velID,bins,aggMeanDBZ,aggStdDBZ,strmPos,
               aggMagThrsh=5,sfcVel=1.5,mskAzBuf_strm=5,mskAzBuf_clr=5,
               mskDifThrsh=10,altAdj=0,zCsumThrsh=10,svExtra=False):
    """
    Function for identifying and masking gates contaminated by the intersection of the
    radar beam with the surface. First, the radar altitude (AGL) is used to identify and 
    mask gates which would reside at or below the surface, assuming level terrain. Next,
    the mean and standard deviation of the reflectivity at each gate (as calculated and
    saved by the `genAgg` method) are added together to produce a clear air masking field.
    Anywhere that field is equal to or greater than some threshold (given by 
    `aggMagThrsh`) becomes a candidate for surface masking. Candidate gates must also
    exist below the radar horizon, with some azimuths above horizon included based on the
    value of `mskAzBuf`. Effects of attenuation are also considered, with gates eligible
    for masking only if the cumulative sum of Z along a given radial is low enough 
    (controlled by the `zCsumThrsh` parameter). Finally, if all above conditions are
    satisfied, the radial velocity associated with a candidate surface gate must be near
    zero (controlled by the `sfcVel` parameter) (the current version of this function will
    also allow masking if the given gate does not have a valid velocity associated with it).
    Enabling the `svExtra` option will output intermediate versions of the reflectivity
    and velocity fields after each of the main steps, as well as difference fields which
    are useful for seeing the exact location and magnitude of any gates masked by these
    operations.
    
    
    Parameters
    ----------
    radar : Py-ART radar object
        Radar object containing the reflectivity and velocity fields to be processed.
    dbzID : string
        Name of the reflectivity variable to be processed.
    velID : string
        Name of the velocity variable to be processed.
    bins : float array
        Array containing azimuth angles used in the generation of the aggregate radar
        fields in `genAgg`. Easiest to generate this array using `np.arange(azPerSwp)`.
        This is used to determine which radials in the current sweep correspond to those
        in the aggregate masking fields.
    aggMeanDBZ : float array
        2D array containing the mean reflectivity at each gate as calculated from the 
        aggregate of all sweeps included in the call to `genAgg`. Should be available 
        within the 'aggCfrad.[...].nc' file.
    aggStdDBZ : float array
        2D array containing the standard deviation of reflectivity at each gate as  
        calculated from the aggregate of all sweeps included in the call to `genAgg`.  
        Should be available within the 'aggCfrad.[...].nc' file.
    strmPos : single character string
        Side of the aircraft heading where the target storm/feature is located. Used in
        application of variable azimuth inclusion/exclusion ranges dependent on where the
        storm is located in radar antenna space. 'L' = storm to the left of aircraft 
        heading, 'R' = storm is to the right of aircraft heading.
    aggMagThrsh : float, optional
        The aggregate masking variable (which is defined as aggMeanDBZ + aggStdDBZ) 
        must be at least this value at a target gate for the gate to be a candidate for
        surface masking. Defaults to 5 dBZ, which seems to work well with at least the TDR 
        data from TORUS19.
    sfcVel : float, optional
        Absolute value of radial velocity (m/s) which a candidate gate for surface masking must 
        be less than or equal to ultimately be masked (i.e., candidate gate with velocity 
        `v` can be masked if this condition is met: `-sfcVel` ≥ `v` ≥ `+sfcVel`).
        Default value is 1.5 m/s, which tends to be a bit conservative. Experimentation
        with this parameter is highly encouraged.
    mskAzBuf_strm : float, optional
        Degrees azimuth above the radar horizon to be searched for surface-contaminated 
        gates on the storm side of the radar. For example, if the storm is to the left of
        the radar, radials with azimuths `azmth` satisfying the following condition will 
        be considered (and the opposite would apply if storm is to the right of the radar):
        (90 - mskAzBuf_clr) ≥ azmth ≥ (270 + mskAzBuf_strm). Defaults to 5 degrees.
    mskAzBuf_clr : float, optional
        Degrees azimuth above the radar horizon to be searched for surface-contaminated 
        gates on the clear air side of the radar. For example, if the storm is to the right of
        the radar, radials with azimuths `azmth` satisfying the following condition will 
        be considered (and the opposite would apply if storm is to the left of the radar):
        (90 - mskAzBuf_strm) ≥ azmth ≥ (270 + mskAzBuf_clr). Defaults to 5 degrees.
    mskDifThrsh : float, optional
        When the difference between reflectivity at a given gate and the reflectivity at 
        the same gate within the clear air aggregate exceeds this threshold value AND the 
        radial velocity is within the range defined by ±sfcVel, the gate will be unmasked 
        if originally flagged by earlier steps. Defaults to 10 dBZ.
    altAdj : float, optional
        Altitude adjustment (meters) to apply to the aircraft altitude for the purposes of
        modifying the geometry-based surface removal. Deafults to 0 (no adjustment).
    zCsumThrsh : float, optional
        The cumulative sum (radial) of Z/(1x10^6) at a target gate must be less than this
        value to be considered for surface masking. Surface returns tend to become less
        significant in regions of greater precipitation (due to attenuation). The default
        value of 10 (which is equivalent to a cumulative sum of Z of 10x10^6 mm^6 m^-3)
        seems to work well with at least the TDR data from TORUS19. If `svExtra` is 
        `True`, the cumulative sum of Z/(1x10^6) will be save to the radar object, which
        can then be used to help identify an appropriate threshold if the default is not
        satisfactory.
    svExtra : bool, optional
        Flag to enable output of intermediate/additional variables useful in assessing
        the performance of this editing step.
    
    
    
    Returns
    -------
    radar : Py-ART radar object
        Py-ART radar object containing original and QC'd data. This radar object can then be used for
        additional analysis or plotting, as well as output to finalized cfRadial files.
    """

    swpAngl = radar.fixed_angle['data'].data[0]

    # Pull in DBZ and VEL and mask regions furthest from radar
    az_raw = radar.rotation['data'].data
    roll = radar.roll['data'].data
    azmth = az_raw + roll
    azmth[azmth < 0] += 360
    azmth[azmth > 360] -= 360
    rng = radar.range['data'].data

    DBZ = radar.fields[dbzID]['data'].copy()
    VEL = radar.fields[velID]['data'].copy()

    
#=================================#
# Apply geometric surface ID/mask #
#=================================#
    geoSfcMsk = np.full_like(DBZ.mask,False)

    radAlt = radar.altitude_agl['data'][0]
    radAlt = radAlt + altAdj

    for aIx,a in enumerate(azmth):
        if a > 90 and a < 270:
            if a < 180:
                sfcRng = radAlt/np.cos(np.deg2rad(180-a))
            #if a > 180:
            #    sfcRng = radAlt/np.cos(np.deg2rad(a-180))
            else:
                sfcRng = radAlt/np.cos(np.deg2rad(180-a))

            if sfcRng <= np.max(rng):
                sfcIx = np.where(rng >= sfcRng)[0][0]
                geoSfcMsk[aIx,sfcIx:] = True


    dbz1_geo = np.ma.masked_array(data=DBZ.data,mask=np.logical_or(geoSfcMsk,DBZ.mask),copy=True)
    #vel1_geo = np.ma.masked_array(data=VEL.data,mask=np.logical_or(geoSfcMsk,VEL.mask),copy=True)
    vel1_geo = np.ma.masked_array(data=VEL.data,mask=np.logical_or(dbz1_geo.mask,VEL.mask),copy=True)
    
    Z_cSum = rayZCsum(dbz1_geo)
        
    
    if svExtra:
        radar.add_field_like(dbzID,dbzID + '_1g',dbz1_geo,replace_existing=True)
        radar.add_field_like(velID,velID + '_1g',vel1_geo,replace_existing=True)
        
        dDBZ_geo = DBZ.filled(fill_value=0)-dbz1_geo.filled(fill_value=0)
        radar.add_field_like(dbzID,dbzID + '_1gDd',dDBZ_geo,replace_existing=True)
        
        radar.add_field_like(dbzID,'Z_cSum',Z_cSum,replace_existing=True)
        

    
#===========================#
# Apply clear air (CA) mask #
#===========================#
    caSfcMsk = np.full_like(DBZ.mask,False)
    loZcSmsk = np.full_like(DBZ.mask,False)
    hiZcSmsk = np.full_like(DBZ.mask,False)
    azMsk = np.full_like(DBZ.mask,False)
    agThrMsk = np.full_like(DBZ.mask,False)
    obExcdMsk = np.full_like(DBZ.mask,False)
    loVelMsk = np.full_like(DBZ.mask,False)
    obMskDifMsk = np.full_like(DBZ.mask,False)
    
    # Experimental masks
    #loloVelMsk = np.full_like(DBZ.mask,False)
    #strmPosAzMsk = np.full_like(DBZ.mask,False)

    # Determine which bins (and the order therein) correspond
    # to the current sweep azimuths
    newIx = (np.digitize(azmth,bins,right=False))-1

    # Generate views of the aggregate mean and std. dev. arrays with
    # only azimuths corresponding to the current sweep
    aggMeanDBZ_v = aggMeanDBZ.copy()[newIx,:]
    aggStdDBZ_v = aggStdDBZ.copy()[newIx,:]

    # Create (and optionally save) the final aggregate clear air masking field
    aggCmb_v = aggMeanDBZ_v + aggStdDBZ_v
    if svExtra:
        radar.add_field_like(dbzID,dbzID + '_1ag',aggCmb_v,replace_existing=True)

    #hiAzIx = np.where(np.logical_and(azmth <= 90, azmth >= 270))[0]
    
    # Find gates where the cumulative Z sum is less than the associated threshold
    loZcSumIx = np.where(Z_cSum < zCsumThrsh)
    loZcSmsk[loZcSumIx] = True
    
    # Find gates where the cumulative Z sum is greater than the associated threshold
    hiZcSumIx = np.where(Z_cSum > zCsumThrsh)
    hiZcSmsk[hiZcSumIx] = True
    
    # Find gates where azimuth is within our targeting range
    if strmPos == 'L':
        azIx = np.where(np.logical_and(azmth >= (90 - mskAzBuf_clr), azmth <= (270 + mskAzBuf_strm)))[0]
    elif strmPos == 'R':
        azIx = np.where(np.logical_and(azmth >= (90 - mskAzBuf_strm), azmth <= (270 + mskAzBuf_clr)))[0]
    azMsk[azIx,:] = True
    
    # Find gates where the aggregate masking field is greater than the
    #   associated threshold
    aggThrIx = np.where(aggCmb_v >= aggMagThrsh)
    agThrMsk[aggThrIx] = True
    
    # Find gates where reflectivity exceeds that of the aggregate masking field
    obExcdIx = np.where(dbz1_geo > aggCmb_v)
    obExcdMsk[obExcdIx] = True
    
    # Find gates where velocity is lower than `sfcVel` or where velocity is masked
    loVelIx = np.where(np.logical_or(np.logical_and(vel1_geo >= -1*sfcVel, vel1_geo <= sfcVel),vel1_geo.mask))
    loVelMsk[loVelIx] = True
    
    
    # Experimental masks
    obMskDifIx = np.where((dbz1_geo - aggCmb_v) > mskDifThrsh)
    obMskDifMsk[obMskDifIx] = True
    
    #loloVelIx = np.where(np.logical_or(np.logical_and(vel1_geo >= -1, vel1_geo <= 1),vel1_geo.mask))
    #loloVelMsk[loloVelIx] = True

    #if strmPos == 'L':
    #    strmPosAzIx = np.where(np.logical_and(azmth >= 180, azmth <= 359))[0]
    #elif strmPos == 'R':
    #    strmPosAzIx = np.where(np.logical_and(azmth >= 0, azmth <= 180))[0]
    #strmPosAzMsk[strmPosAzIx] = True
    
    
    # Original masking set
    #caSfcMsk[np.logical_and(np.logical_and(azMsk,np.logical_and(agThrMsk,loVelMsk)),loZcSmsk)] = True
    #caSfcMsk[np.logical_and(obExcdMsk,np.logical_not(loVelMsk))] = False
    
    # New (11/12/19) masking set
    caSfcMsk[np.logical_and(np.logical_and(azMsk,np.logical_and(agThrMsk,loVelMsk)),loZcSmsk)] = True # Mask gates near/below radar horizon if masking field exceeds threshold, velocity is within surface range, and cumulative Z sum is lower than the threshold
    caSfcMsk[np.logical_and(np.logical_and(hiZcSmsk,loVelMsk),azMsk)] = True # Mask gates with high cumulative Z sum and velocity within surface range IF near/below radar horizon
    caSfcMsk[np.logical_and(obExcdMsk,np.logical_not(loVelMsk))] = False # Unmask any gates with higher DBZ than same gate in mask IF velocity is greater than surface range
    caSfcMsk[np.logical_and(obMskDifMsk,loVelMsk)] = False # Unmask any gates with higher DBZ (by at least mskDifThrsh) than same gate in mask IF velocity is within surface range
        

    dbz2_m = np.ma.masked_array(data=DBZ.data,mask=np.logical_or(caSfcMsk,dbz1_geo.mask),copy=True)
    vel2_m = np.ma.masked_array(data=VEL.data,mask=np.logical_or(dbz2_m.mask,vel1_geo.mask),copy=True)
    radar.add_field_like(dbzID,dbzID,dbz2_m,replace_existing=True)
    radar.add_field_like(velID,velID,vel2_m,replace_existing=True)

    if svExtra:
        dDBZ_m = dbz1_geo.filled(fill_value=0)-dbz2_m.filled(fill_value=0)
        radar.add_field_like(dbzID,dbzID + '_1mDd',dDBZ_m,replace_existing=True)

    return radar


def radFilt(radar,thrshVarID,dbzID,velID,thrshFillVal,filtStdX=5,filtStdY=1,smthThrsh=0,
            retainShrtPulse=False,svExtra=False):
    """
    
    
    
    Parameters
    ----------
    radar : Py-ART radar object
        Radar object containing the reflectivity and velocity fields to be processed.
    thrshVarID : string
        Name of variable to be used as the thresholding field. This variable will be
        smoothed and then compared against `thrshFillVal` to threshold reflectivity
        and velocity.
    dbzID : string
        Name of the reflectivity variable to be processed.
    velID : string
        Name of the velocity variable to be processed.
    thrshFillVal : float
        Value to fill any missing/NaN gates with within thresholding field as a part of 
        the smoothing process. Typically, this value should be among the lowest expected
        values for the chosen thresholding field (e.g., -10 works well when reflectivity
        is used as the thresholding field).
    filtStdX : float, optional
       Standard deviation of the Gaussian kernel in X used to create the smoothed masking field.
       The smaller this value is, the sharper the filter response. A larger value of filtStdX
       (i.e., > 5) and smaller value of filtStdY (i.e., < 2) seems to work well on reflectivity.
    filtStdY : float, optional
       Standard deviation of the Gaussian kernel in Y used to create the smoothed masking field.
       The smaller this value is, the sharper the filter response. A larger value of filtStdX
       (i.e., ≥ 5) and smaller value of filtStdY (i.e., < 2) seems to work well on reflectivity.
    smthThrsh : float, optional
        Any gate where the smoothed reflectivity field is less than `smthThrsh` (defaults 
        to 0) is masked in the final output fields.
    retainShrtPulse : bool, optional
        Flag to control filtering of reflectivity for gates obtained using the short pulse  
        (of a dual-PRT pair). This specifically applies to TORUS TDR data between 450 m  
        and ~975 m range. If similar behavior is desired for other radars/ranges, this 
        function will need to be modified accordingly. Defaults to False.
    svExtra : bool, optional
        Flag to enable output of intermediate/additional variables useful in assessing
        the performance of this editing step.
    
    
    
    Returns
    -------
    radar : Py-ART radar object
        Py-ART radar object containing original and QC'd data. This radar object can then be used for
        additional analysis or plotting, as well as output to finalized cfRadial files.
    """
    DBZ = radar.fields[dbzID]['data'].copy()
    VEL = radar.fields[velID]['data'].copy()
    
    if thrshVarID == dbzID:
        thrshVar = DBZ.copy()
    elif thrshVarID == velID:
        thrshVar = VEL.copy()
    else:
        thrshVar = radar.fields[thrshVarID]['data'].copy()
    
    g2dKernel = Gaussian2DKernel(x_stddev=filtStdX,y_stddev=filtStdY)
    thrsh_smth = convolve(thrshVar,g2dKernel,fill_value=thrshFillVal,nan_treatment='fill',
                        boundary='extend',preserve_nan=True)
    thrsh_smth = np.ma.masked_invalid(thrsh_smth)
        

    smthFltMsk = np.full_like(DBZ.mask,False)
    smthFlgIx = np.ma.where(thrsh_smth < smthThrsh)
    smthFltMsk[smthFlgIx] = True

    dbz_filt = np.ma.masked_array(data=DBZ.data,mask=np.logical_or(smthFltMsk,DBZ.mask),copy=True)
    vel_filt = np.ma.masked_array(data=VEL.data,mask=np.logical_or(dbz_filt.mask,VEL.mask),copy=True)
    
    if retainShrtPulse:
        dbz_filt[:,6:14] = DBZ[:,6:14]
        vel_filt[:,6:14] = VEL[:,6:14]
    
    radar.add_field_like(dbzID,dbzID,dbz_filt,replace_existing=True)
    radar.add_field_like(velID,velID,vel_filt,replace_existing=True)

    
    if svExtra:
        radar.add_field_like(thrshVarID,thrshVarID + '_sm',thrsh_smth,replace_existing=True)#.filled(fill_value=-10)
        dDBZ_filt = DBZ.filled(fill_value=0)-dbz_filt.filled(fill_value=0)
        radar.add_field_like(dbzID,dbzID + '_fDd',dDBZ_filt,replace_existing=True)
    
    return radar
    

def caVelChk(radar,dbzID,velID,simVelID,refThrsh=10,velDifThrsh=10,svExtra=False):
    """
    
    
    
    Parameters
    ----------
    radar : Py-ART radar object
        Radar object containing the reflectivity and velocity fields to be processed.
    dbzID : string
        Name of the reflectivity variable to be used for reflectivity thresholding.
    velID : string
        Name of the observed radial velocity variable to be compared with the simulated
        radial velocity field.
    simVelID : string
        Name of the simulated radial velocity field. This field can be generated using the
        Py-ART `HorizontalWindProfile.from_u_and_v` the output of which is then fed to 
        another Py-ART function, `simulated_vel_from_profile`.
    refThrsh : float, optional
        Max value of reflectivity to be considered clear air. Observed velocities will 
        only be compared to simulated radial velocities when the gate reflectivity is
        less than or equal to this value. Defaults to 10 dBZ.
    velDifThrsh : float, optional
        Absolute difference allowed between the observed radial velocity and the 
        simulated radial velocity. When the difference exceeds this value and has a 
        reflectivity less than or equal to `refThrsh`, that gate will be masked in `velID`
        only (reflectivity variable `dbzID` is not modified by this function). Defaults
        to 10 m/s.
    svExtra : bool, optional
        Flag to enable output of intermediate/additional variables useful in assessing
        the performance of this editing step.
    
    
    
    Returns
    -------
    radar : Py-ART radar object
        Py-ART radar object containing original and QC'd data. This radar object can then be used for
        additional analysis or plotting, as well as output to finalized cfRadial files.
    """
    DBZ = radar.fields[dbzID]['data'].copy()
    VEL = radar.fields[velID]['data'].copy()
    try:
        velSim = radar.fields[simVelID]['data'].copy()
    except:
        print('{} does not exist in the radar file. Please check that you have properly entered the simVelID. Exiting...'.format(simVelID))
        
        
    # Get absolute difference between observed and simulated winds
    velDiff = VEL - velSim
    velDiffAbs = np.abs(velDiff)
    
    lowRefMsk = np.full_like(DBZ.mask,False)
    lrFlgIx = np.ma.where(DBZ <= refThrsh)
    lowRefMsk[lrFlgIx] = True
    
    caVelMsk = np.full_like(DBZ.mask,False)
    cavFlgIx = np.ma.where(velDiffAbs > velDifThrsh)
    caVelMsk[cavFlgIx] = True

    clrAdjMsk = np.full_like(DBZ.mask,False)
    clrAdjMsk[np.logical_and(lowRefMsk,caVelMsk)] = True

    vel_clrAdj = np.ma.masked_array(data=VEL.data,mask=np.logical_or(clrAdjMsk,VEL.mask),copy=True)
    dbz_clrAdj = np.ma.masked_array(data=DBZ.data,mask=np.logical_or(clrAdjMsk,DBZ.mask),copy=True)
    
    radar.add_field_like(velID,velID,vel_clrAdj,replace_existing=True)
    radar.add_field_like(dbzID,dbzID,dbz_clrAdj,replace_existing=True)
    
    if svExtra:
        radar.add_field_like(velID,velID + '_simDif',velDiff,replace_existing=True)
        dVEL_filt = VEL.filled(fill_value=0)-vel_clrAdj.filled(fill_value=0)
        radar.add_field_like(velID,velID + '_simDd',dVEL_filt,replace_existing=True)
        
    return radar