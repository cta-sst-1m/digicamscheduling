import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord, AltAz
import numpy as np
import matplotlib.pyplot as plt


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


def intensity(date, location, ra, dec):

    alt = compute_source_position(date, location, ra, dec).alt

    return compute_source_intensity(alt)