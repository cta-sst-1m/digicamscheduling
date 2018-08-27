import astropy.units as u
from scipy.interpolate import interp1d
import numpy as np


def compute_forbidden_zone(alt, az):

    return alt > 45 * u.deg


def interpolate_environmental_limits(alt, az, **kwargs):
    fov = 4.5 * u.deg
    alt = alt + fov
    f = interp1d(az, alt, kind='cubic', **kwargs)

    return f


def is_above_environmental_limits(alt, az, environmental_limits):

    az = az.to(u.deg).value
    alt = alt.to(u.deg).value

    return alt > (environmental_limits(az))


def compute_observability(sun_elevation, moon_alt, moon_ph, use_moon=True):

    moon_elevation = moon_alt.copy()
    moon_phase = moon_ph.copy()

    moon_phase[(moon_elevation < 0 * u.deg)] = 0.
    moon_elevation[(moon_elevation < 0 * u.deg)] = 0 * u.deg
    observability = (sun_elevation < -12 * u.deg)

    if use_moon:

        observability = observability * np.cos(moon_elevation)
        observability = observability * (1. - moon_phase)

    return observability
