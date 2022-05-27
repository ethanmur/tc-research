"""
pyart.util.simulated_vel
========================

Function for creating simulated velocity fields.

.. autosummary::
    :toctree: generated/

    simulated_vel_from_profile

"""

import numpy as np
from scipy.interpolate import interp1d

from ..config import get_metadata, get_field_name
from ..core import antenna_vectors_to_cartesian


def simulated_vel_from_profile(
        radar, profile, interp_kind='linear', sim_vel_field=None, airTail=False):
    """
    Create simulated radial velocities from a profile of horizontal winds.

    Parameters
    ----------
    radar : Radar
        Radar instance which provides the scanning parameters for the
        simulated radial velocities.
    profile : HorizontalWindProfile
        Profile of horizontal winds.
    interp_kind : str, optional
        Specifies the kind of interpolation used to determine the winds at a
        given height. Must be one of 'linear', 'nearest', 'zero', 'slinear',
        'quadratic', or 'cubic'. The the documentation for the SciPy
        scipy.interpolate.interp1d function for descriptions.
    sim_vel_field : str, optional
        Name to use for the simulated velocity field metadata. None will use
        the default field name from the Py-ART configuration file.
    airTail : bool, optional
        Indicate whether the input radar is an airborne tail platform (e.g., NOAA TDR)
        to appropriately adjust the coordinate transformation.

    Returns
    -------
    sim_vel : dict
        Dictionary containing a radar field of simulated radial velocities.

    """
    # parse parameters
    if sim_vel_field is None:
        sim_vel_field = get_field_name('simulated_velocity')


    # radar parameters
    if airTail:
        az,el = antenna_to_track_relative_dpj(radar)
        azimuths = np.deg2rad(az).reshape(-1, 1)
        elevations = np.deg2rad(el).reshape(-1, 1)
        
        xg,yg,zg = antenna_vectors_to_cartesian(radar.range['data'],az,el)
        
        gate_altitudes = np.mean(radar.altitude['data']) + zg
    else:
        azimuths = np.deg2rad(radar.azimuth['data']).reshape(-1, 1)
        elevations = np.deg2rad(radar.elevation['data']).reshape(-1, 1)
        gate_altitudes = radar.gate_altitude['data']

    if isinstance(gate_altitudes, np.ma.MaskedArray):
        gate_altitudes = gate_altitudes.filled(np.nan)

    # prepare wind profile for interpolation
    if isinstance(profile.height, np.ma.MaskedArray):
        height = profile.height.filled(np.nan)
    else:
        height = profile.height

    height_is_not_nan = ~np.isnan(height)
    winds = np.empty((2, len(height)), dtype=np.float64)
    if isinstance(profile.u_wind, np.ma.MaskedArray):
        winds[0] = profile.u_wind.filled(np.nan)
    else:
        winds[0] = profile.u_wind

    if isinstance(profile.v_wind, np.ma.MaskedArray):
        winds[1] = profile.v_wind.filled(np.nan)
    else:
        winds[1] = profile.v_wind

    wind_is_not_nan = np.logical_and(~np.isnan(winds[0]), ~np.isnan(winds[1]))
    no_nans = np.logical_and(height_is_not_nan, wind_is_not_nan)
    height = height[no_nans]
    winds[0] = winds[0][no_nans]
    winds[1] = winds[1][no_nans]
    wind_interp = interp1d(
        height, winds, kind=interp_kind, bounds_error=False)

    # interpolated wind speeds at all gates altitudes
    gate_winds = wind_interp(gate_altitudes)
    gate_u = np.ma.masked_invalid(gate_winds[0])
    gate_v = np.ma.masked_invalid(gate_winds[1])

    # calculate the radial velocity for all gates
    radial_vel = (gate_u * np.sin(azimuths) * np.cos(elevations) +
                  gate_v * np.cos(azimuths) * np.cos(elevations))

    sim_vel = get_metadata(sim_vel_field)
    sim_vel['data'] = radial_vel
    return sim_vel


def antenna_to_track_relative_dpj(radar):
    rotD_raw = radar.rotation['data']
    tiltD = radar.tilt['data']
    pitchD = radar.pitch['data']
    rollD = radar.roll['data']
    driftD = radar.drift['data']
    trackD = radar.heading['data']

    rotD = rotD_raw + rollD
    rotD[rotD < 0] += 360
    rotD[rotD > 360] -= 360

    rot = np.deg2rad(rotD)
    tilt = np.deg2rad(tiltD)
    pitch = np.deg2rad(pitchD)
    drift = np.deg2rad(driftD)
    
    # Direction cosine for x (distance normal to the track)
    x = ((np.cos(rot)*np.sin(drift)*np.sin(pitch)*np.cos(tilt))
         + (np.cos(drift)*np.sin(rot)*np.cos(tilt))
         - (np.sin(drift)*np.cos(pitch)*np.sin(tilt)))

    # Direction cosine for y (distance along track due to fore/aft pointing)
    y = ((-1.0*np.cos(rot)*np.cos(drift)*np.sin(pitch)*np.cos(tilt)) 
         + (np.sin(drift)*np.sin(rot)*np.cos(tilt)) 
         + (np.cos(drift)*np.cos(pitch)*np.sin(tilt)))

    # Direction cosine for z (height)
    z = ((np.cos(pitch)*np.cos(rot)*np.cos(tilt)) 
         + (np.sin(pitch)*np.sin(tilt)))

    # Azimuth (x,y) angle (compass direction) relative to the track
    az_tr = np.rad2deg(np.arctan2(x,y))

    # Azimuth angle relative to North
    azimuths = np.mod((az_tr + trackD), 360.0)

    # Elevation angle from the horizontal (+ or - 90 deg)
    elevations = np.rad2deg(np.arcsin(z))
    
    return azimuths,elevations