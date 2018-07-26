import astropy.units as u
from scipy.interpolate import interp1d


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