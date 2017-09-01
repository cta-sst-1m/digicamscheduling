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

if __name__ == '__main__':

    from digicam_scheduling.core import moon

    sources = [{'ra': 83.633212 * u.deg, 'dec': 22.01446 * u.deg, 'name': 'Crab'}, {'dec': 38.208801 * u.deg, 'name': 'Mrk 421', 'ra': 166.113808 * u.deg, 'flux': 0.3}]
    coordinates_krakow = {'lat': 50.090815 * u.deg, 'lon': 19.887937 * u.deg, 'height': 214.034 * u.m}
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 1, num=1000) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-08-31 12:00')
    source_intensity = np.zeros((len(sources), time_bins.shape[0]))
    date = start_date + time_bins

    for i in range(len(sources)):

        altaz = compute_source_position(date=date + time_bins, location=location, ra=sources[i]['ra'], dec=sources[i]['dec'])



    plt.show()