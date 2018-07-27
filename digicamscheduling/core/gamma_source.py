import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
import numpy as np
from digicamscheduling.utils.decorator import timeit


# @timeit
def compute_source_position(date, location, ra, dec):

    altaz_coordinates = AltAz(obstime=date, location=location)
    source_location = SkyCoord(ra=ra, dec=dec, frame='icrs')
    source_altaz = source_location.transform_to(altaz_coordinates)

    return source_altaz


def compute_source_intensity(alt, trees=True):

    if trees:

        return np.sin(alt) * (alt > 45 * u.deg) * (alt < 85 * u.deg)

    else:

        return np.sin(alt) * (alt > 0 * u.deg)


def intensity(date, location, ra, dec, altaz_moon=None):

    altaz_source = compute_source_position(date, location, ra, dec)

    if altaz_moon is None:

        return compute_source_intensity(altaz_source.alt)

    else:

        moon_distance = compute_moon_distance(altaz_source=altaz_source, altaz_moon=altaz_moon)
        return compute_source_intensity(altaz_source.alt) * (moon_distance > 15 * u.deg)


def compute_moon_distance(altaz_source, altaz_moon):

    return altaz_source.separation(altaz_moon)



def compute_observatility(sun_alt, moon_alt, moon_phase):
    moon_elevation[moon_elevation < 0 * u.deg] = 0
    observability = (sun_elevation < -12 * u.deg)
    observability = observability * np.cos(moon_elevation)
    observability *= (
    1 - moon_phase * (moon_elevation > 0 * u.deg).astype(float))