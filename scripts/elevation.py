import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source
import matplotlib.pyplot as plt


if __name__ == '__main__':

    sources_filename = '../digicamscheduling/config/' + 'fact_catalog.txt'
    location_filename = '../digicamscheduling/config/' + 'location_krakow.txt'

    sources = reader.read_catalog(sources_filename)
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 0.5, num=12*2 + 1) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-08-31 20:00')

    source_elevation = np.zeros((len(sources), time_bins.shape[0])) * u.deg
    source_azimut = np.zeros((len(sources), time_bins.shape[0])) * u.deg

    for i, time in enumerate(time_bins):

        date = start_date + time

        for j, source in enumerate(sources):

            temp = gamma_source.compute_source_position(date=date, location=location, ra=source['ra'], dec=source['dec'])
            source_elevation[j, i] = temp.alt
            source_azimut[j, i] = temp.az

    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()

    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')

    for j, source in enumerate(sources):

        t = (start_date + time_bins).plot_date
        alt = source_elevation[j]
        az = source_azimut[j]
        axis_1.plot_date(t, alt, label=sources[j]['name'], linestyle='-', marker='None')
        axis_2.plot_date(t, az, label=sources[j]['name'], linestyle='-', marker='None')
        axis_3.plot(az.to('radian'), (90 * u.deg - alt), label=sources[j]['name'])

    axis_1.set_xlabel('UTC time')
    axis_2.set_xlabel('UTC time')
    axis_1.set_ylabel('elevation [deg]')
    axis_2.set_ylabel('azimuth [deg]')
    axis_1.set_ylim([0, 90])
    axis_3.set_rmax(90)
    axis_1.legend(loc='best')
    axis_2.legend(loc='best')
    axis_3.legend(loc='best')
    plt.gcf().autofmt_xdate()
    plt.show()


