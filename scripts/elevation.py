import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source
import digicamscheduling.display.plot as display
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
    date = start_date + time_bins

    source_elevation = np.zeros((len(sources), time_bins.shape[0])) * u.deg
    source_azimut = np.zeros((len(sources), time_bins.shape[0])) * u.deg

    for i, source in enumerate(sources):

        temp = gamma_source.compute_source_position(date=date, location=location, ra=source['ra'], dec=source['dec'])
        source_elevation[i] = temp.alt
        source_azimut[i] = temp.az

    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()

    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')

    for j, source in enumerate(sources):

        alt = source_elevation[j]
        az = source_azimut[j]
        display.plot_azimuth(date, az, axis=axis_1, label=source['name'])
        display.plot_elevation(date, alt, axis=axis_2, label=source['name'])
        display.plot_trajectory(az, alt, axis=axis_3, label=source['name'])


    plt.show()


