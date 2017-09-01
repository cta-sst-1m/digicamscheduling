from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import core.sun as sun
import core.moon as moon
import core.gamma_source as gamma_source

if __name__ == '__main__':

    coordinates_krakow = {'lat': 50.090763 * u.deg, 'lon': 19.887956 * u.deg, 'height': 230 * u.m}
    location = EarthLocation(**coordinates_krakow)
    krakow_utc = 2. * u.hour

    sky_coordinates = SkyCoord.from_name('crab')
    ra = sky_coordinates.ra
    dec = sky_coordinates.dec

    time_intervals = np.linspace(0, 7, num=2000) * u.day
    crab_intensity = np.zeros(time_intervals.shape[0])
    moon_intensity = np.zeros(time_intervals.shape[0])
    sun_intensity = np.zeros(time_intervals.shape[0])

    time = Time('2017-08-27 00:00')
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