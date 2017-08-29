from astropy.coordinates import AltAz, get_sun
from astropy import units as u
# import numpy as np


def compute_sun_position(date, location):

    altaz_coordinates = AltAz(obstime=date, location=location)
    sun_altaz = get_sun(date).transform_to(altaz_coordinates)

    return sun_altaz


def compute_sun_intensity(alt):

    return alt >= - 18 * u.deg


def intensity(date, location):

    alt = compute_sun_position(date, location).alt

    return compute_sun_intensity(alt)