from astropy.coordinates import AltAz, get_sun, EarthLocation
from astropy.time import Time
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt


def compute_sun_position(date, location):

    altaz_coordinates = AltAz(obstime=date, location=location)
    sun_altaz = get_sun(date).transform_to(altaz_coordinates)

    return sun_altaz


def compute_sun_intensity(alt):

    return alt >= - 18 * u.deg


def compute_night(alt, type='astronomical'):

    night = (alt <= -12 * u.deg)

    return night


def intensity(date, location):

    alt = compute_sun_position(date, location).alt

    return compute_sun_intensity(alt)


if __name__ == '__main__':

    coordinates_krakow = {'lat': 50.090763 * u.deg, 'lon': 19.887956 * u.deg,
                          'height': 230 * u.m}
    location = EarthLocation(**coordinates_krakow)

    time_intervals = np.linspace(0, 7, num=500) * u.day
    alt = np.zeros(time_intervals.shape[0]) * u.deg
    az = np.zeros(time_intervals.shape[0]) * u.deg
    intensity_sun = np.zeros(time_intervals.shape[0])

    time = Time('2017-08-27 00:00')
    date_for_plot = time + time_intervals

    for i in range(time_intervals.shape[0]):
        date = time + time_intervals[i]
        temp = compute_sun_position(date, location=location)
        alt[i] = temp.alt.to('deg')
        az[i] = temp.az.to('deg')
        intensity_sun[i] = intensity(date=date, location=location)

    # sky_darkness = (1 - compute_moon_intensity(alt, phases)) *
    # (1 - compute_sun_intensity(alt_sun))

    fig, axis_1 = plt.subplots()
    axis_1.plot_date(date_for_plot.plot_date, intensity_sun,
                     label='sun intensity', color='k', linestyle='--',
                     marker='None')
    # axis_1.plot_date(date_for_plot.plot_date, sky_darkness,
    # label='sky darkness', color='k', linestyle='--', marker='None')
    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[a.u.]', color='k')
    axis_1.set_ylim([0, 1.2])
    axis_1.tick_params('y', colors='k')
    plt.legend(loc='upper left')
    plt.gcf().autofmt_xdate()
    axis_2 = axis_1.twinx()
    axis_2.plot_date(date_for_plot.plot_date, alt,
                     color='r', label='Sun', linestyle='-', marker='None')
    # axis_2.plot_date(date_for_plot.plot_date, alt_sun,
    #  color='r', label='Sun', linestyle='-', marker='None' )
    axis_2.set_ylabel('altitude [deg]', color='r')
    axis_2.set_ylim([0, 90])
    axis_2.tick_params('y', colors='r')
    plt.legend(loc='upper right')
    plt.show()


