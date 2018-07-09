import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy import pyasl
from astropy.coordinates import AltAz, EarthLocation, get_moon
from astropy.time import Time


def compute_moon_phase(date):

    julian_date = date.jd

    return pyasl.moonphase(julian_date)


def compute_moon_position(date, location):

    altaz_coordinates = AltAz(obstime=date, location=location)
    moon_altaz = get_moon(date).transform_to(altaz_coordinates)

    return moon_altaz


def compute_moon_intensity(alt, phase):

    return phase * np.sin(alt) * (alt > 0. * u.deg)


def intensity(date, location):

    alt = compute_moon_position(date, location).alt
    phase = compute_moon_phase(date)

    return compute_moon_intensity(alt, phase)


if __name__ == '__main__':

    import sun

    coordinates_krakow = {'lat': 50.090763 * u.deg, 'lon': 19.887956 * u.deg, 'height': 230 * u.m}
    location = EarthLocation(**coordinates_krakow)
    krakow_utc = 2. * u.hour

    time_intervals = np.linspace(0, 7, num=500) * u.day
    phases = np.zeros(time_intervals.shape[0])
    alt = np.zeros(time_intervals.shape[0]) * u.deg
    az = np.zeros(time_intervals.shape[0]) * u.deg
    intensity_moon = np.zeros(time_intervals.shape[0])
    intensity_sun = np.zeros(time_intervals.shape[0])
    sky_darkness = np.zeros(time_intervals.shape[0])


    time = Time('2017-08-27 00:00')
    date_for_plot = time + time_intervals

    for i in range(time_intervals.shape[0]):
        date = time + time_intervals[i]
        phases[i] = compute_moon_phase(date)
        temp = compute_moon_position(date, location=location)
        alt[i] = temp.alt.to('deg')
        az[i] = temp.az.to('deg')
        intensity_moon[i] = intensity(date=date, location=location)
        intensity_sun[i] = sun.intensity(date=date, location=location)

    fig, axis_1 = plt.subplots()
    axis_1.plot_date(date_for_plot.plot_date, phases, label='moon fraction', color='k', linestyle='-', marker='None')
    axis_1.plot_date(date_for_plot.plot_date, intensity_moon, label='moon intensity', color='k', linestyle='--', marker='None')
    # axis_1.plot_date(date_for_plot.plot_date, sky_darkness, label='sky darkness', color='k', linestyle='--', marker='None')
    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[a.u.]', color='k')
    axis_1.set_ylim([0, 1.2])
    axis_1.tick_params('y', colors='k')
    plt.legend(loc='upper left')
    plt.gcf().autofmt_xdate()
    axis_2 = axis_1.twinx()
    axis_2.plot_date(date_for_plot.plot_date, alt, color='r', label='Moon', linestyle='-', marker='None')
    # axis_2.plot_date(date_for_plot.plot_date, alt_sun, color='r', label='Sun', linestyle='-', marker='None' )
    axis_2.set_ylabel('altitude [deg]', color='r')
    axis_2.set_ylim([0, 90])
    axis_2.tick_params('y', colors='r')
    plt.legend(loc='upper right')

    """
    axis_3 = axis_1.twinx()
    axis_3.plot(time_intervals, az, color='g')
    axis_3.set_ylabel('azimuth [deg]', color='g')
    axis_3.tick_params('y', colors='g')
    """


    sky_darkness = (1 - intensity_moon) * (1 - intensity_sun)

    fig, axis_1 = plt.subplots()
    axis_1.plot_date(date_for_plot.plot_date, sky_darkness, label='sky darkness', color='k', linestyle='-', marker='None')
    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[a.u.]', color='k')
    axis_1.set_ylim([0, 1.2])
    axis_1.tick_params('y', colors='k')
    plt.gcf().autofmt_xdate()
    plt.legend(loc='best')

    print(np.max(alt), np.min(alt))
    print(np.max(az), np.min(az))

    mask =(alt <=90 * u.deg) * (alt >= 0 * u.deg)

    plt.figure()
    axis_polar = plt.subplot(111, projection='polar')
    axis_polar.plot(az[mask], alt[mask])
    axis_polar.set_rmax(90)
    plt.show()
