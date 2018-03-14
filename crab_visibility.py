from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from digicamscheduling.core import sun, moon, gamma_source
from digicamscheduling.io import reader


def crab_visibility(location_filename='digicamscheduling/config/location_krakow.txt'):

    # coordinates_krakow = {'lat': 50.090763 * u.deg, 'lon': 19.887956 * u.deg, 'height': 230 * u.m}
    coordinates_krakow = reader.read_location(location_filename)
    location = EarthLocation(**coordinates_krakow)

    sky_coordinates = SkyCoord.from_name('crab')
    ra = sky_coordinates.ra
    dec = sky_coordinates.dec

    time_intervals = np.linspace(0, 365, num=365) * u.day
    crab_intensity = np.zeros(time_intervals.shape[0])
    moon_intensity = np.zeros(time_intervals.shape[0])
    sun_intensity = np.zeros(time_intervals.shape[0])

    time = Time('2018-01-01 00:00')
    date_for_plot = time + time_intervals

    for i in range(time_intervals.shape[0]):

        date = time + time_intervals[i]

        sun_intensity[i] = sun.intensity(date, location)
        moon_intensity[i] = moon.intensity(date, location)
        crab_intensity[i] = gamma_source.intensity(date, location, ra=ra, dec=dec)

    crab_visibility = (1 - sun_intensity) * (1 - moon_intensity) * crab_intensity

    fig, axis_1 = plt.subplots()
    axis_1.plot_date(date_for_plot.plot_date, crab_visibility, label='Crab', color='k', linestyle='-', marker='None')
    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[u.a]')
    axis_1.set_ylim([0, 1.2])
    plt.legend(loc='best')
    plt.gcf().autofmt_xdate()
    plt.show()


if __name__ == '__main__':

    crab_visibility()