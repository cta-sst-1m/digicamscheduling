import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
import digicamscheduling.core.scheduler as scheduler
import digicamscheduling.core.sun as sun
import digicamscheduling.core.moon as moon
import digicamscheduling.core.gamma_source as gamma_source
from digicamscheduling.io import writer
from digicamscheduling.io import reader

import matplotlib.pyplot as plt


if __name__ == '__main__':

    output_directory = '../digicamscheduling/samples/'
    units_output = 'deg'
    catalog_filename = '../digicamscheduling/config/fact_catalog.txt'
    location_filename = '../digicamscheduling/config/location_krakow.txt'

    sources = reader.read_catalog(filename=catalog_filename)
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 8, num=8*24*2 + 1) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-09-11 12:00')
    date = start_date + time_bins

    sun_intensity = sun.intensity(date=date, location=location)
    moon_intensity = moon.intensity(date=date, location=location)
    source_intensity = np.zeros((len(sources), time_bins.shape[0]))

    for i, source in enumerate(sources):

        source_intensity[i] = gamma_source.intensity(date=date, location=location, ra=source['ra'], dec=source['dec'])

    sources_visibility = scheduler.compute_source_visibility(source_intensity, sun_intensity, moon_intensity)

    availability, schedule_quality = scheduler.find_quality_schedule(sources_visibility=sources_visibility)

    writer.write_schedule(schedule=schedule_quality, sources=sources, start_date=start_date, filename=output_directory + 'quality_schedule_%s.txt' % units_output, time_bins=time_bins, units=units_output)

    fig_1, axis_1 = plt.subplots()
    axis_1.set_title('Quality schedule')
    fig_2, axis_2 = plt.subplots()
    axis_2.set_title('Priority schedule')
    fig_3, axis_3 = plt.subplots()
    axis_3.set_title('Dynamic quality-priority schedule')
    for i, source in enumerate(sources):

        if np.any(schedule_quality[i] > 0):
            axis_1.plot_date(date.plot_date, schedule_quality[i], label=source['name'], linestyle='-', marker='None')

    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[occupancy]')
    axis_1.set_ylim([0, 1.2])
    axis_1.legend(loc='best')
    plt.gcf().autofmt_xdate()
    plt.show()
