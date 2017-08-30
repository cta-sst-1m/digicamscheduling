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


def compute_source_intensity(alt):

    return np.sin(alt) * (alt > 0 * u.deg)


def intensity(date, location, ra, dec):

    alt = compute_source_position(date, location, ra, dec).alt

    return compute_source_intensity(alt)

if __name__ == '__main__':

    sources = [{'ra': 11.074266 * u.deg, 'dec': 38.208801 * u.deg, 'name': 'Mrk 421'},
               {'ra': 5.575539 * u.deg, 'dec': 22.014500 * u.deg, 'name': 'Crab'},
               {'ra': 16.897867 * u.deg, 'dec': 39.760201 * u.deg, 'name': 'Mrk 501'}]

    coordinates_krakow = {'lat': 50.090815 * u.deg, 'lon': 19.887937 * u.deg, 'height': 214.034 * u.m}
    location = EarthLocation(**coordinates_krakow)


    time_bins = np.linspace(0, 7, num=500) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-08-27 00:00')
    source_intensity = np.zeros((len(sources), time_bins.shape[0]))
    date_for_plot = start_date + time_bins


    for i, time in enumerate(time_bins):

        date = start_date + time

        for j, source in enumerate(sources):

            source_intensity[j, i] = intensity(date=date, location=location, ra=source['ra'], dec=source['dec'])

    fig, axis_1 = plt.subplots()

    for i in range(len(sources)):

        axis_1.plot_date(date_for_plot.plot_date, source_intensity[i], label='%s intensity' % (sources[i]['name']), linestyle='-', marker='None')
        axis_1.set_xlabel('UTC time')
        axis_1.set_ylabel('[a.u.]', color='k')
        axis_1.set_ylim([0, 1.2])
        axis_1.tick_params('y', colors='k')
        plt.legend(loc='upper left')
        plt.gcf().autofmt_xdate()

    plt.show()