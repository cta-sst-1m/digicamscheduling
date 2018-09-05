from PyAstronomy import pyasl
from astropy.coordinates import AltAz, get_moon


def compute_moon_phase(date):

    julian_date = date.jd

    return pyasl.moonphase(julian_date)


def compute_moon_position(date, location):

    altaz_coordinates = AltAz(obstime=date, location=location)
    moon_altaz = get_moon(date).transform_to(altaz_coordinates)

    return moon_altaz
