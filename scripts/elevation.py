import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from core.scheduler import compute_source_visibility, find_priority_schedule, compute_schedule_efficiency, find_quality_schedule, find_dynamic_priority_quality_schedule
import core.sun as sun
import core.moon as moon
import core.gamma_source as gamma_source
import matplotlib.pyplot as plt
from inout.reader import read_catalog, read_location


if __name__ == '__main__':

    sources_filename = 'config/' + 'fact_catalog.txt'
    location_filename = 'config/' + 'location_krakow.txt'

    sources = read_catalog(sources_filename)
    coordinates_krakow = read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 0.5, num=12*2 + 1) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-08-31 20:00')

    source_elevation = np.zeros((len(sources), time_bins.shape[0])) * u.deg
    source_azimut = np.zeros((len(sources), time_bins.shape[0])) * u.deg

    for i, time in enumerate(time_bins):

        date = start_date + time

        for j, source in enumerate(sources):

            if source['name'] == '1ES 2344+51.4':

                temp = gamma_source.compute_source_position(date=date, location=location, ra=source['ra'], dec=source['dec'])
                source_elevation[j, i] = temp.alt
                source_azimut[j, i] = temp.az

    fig_1 = plt.figure()
    fig_2 = plt.figure()

    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)

    for j, source in enumerate(sources):


        if source['name'] == '1ES 2344+51.4':
            x = (start_date + time_bins).plot_date
            y = source_elevation[j]
            z = source_azimut[j]

            axis_1.plot_date(x, y, label=sources[j]['name'], linestyle='-', marker='None')
            axis_2.plot_date(x, z, label=sources[j]['name'], linestyle='-', marker='None')

    axis_1.set_xlabel('UTC time')
    axis_2.set_xlabel('UTC time')
    axis_1.set_ylabel('elevation [deg]')
    axis_2.set_ylabel('azimuth [deg]')
    # axis_1.set_ylim([0, 1.1])
    plt.legend(loc='best')
    plt.gcf().autofmt_xdate()
    plt.show()