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
from digicamscheduling.utils import time

import matplotlib.pyplot as plt


if __name__ == '__main__':

    output_directory = '../digicamscheduling/samples/'
    units_output = 'deg'
    catalog_filename = '../digicamscheduling/config/fact_catalog.txt'
    location_filename = '../digicamscheduling/config/location_krakow.txt'

    sources = reader.read_catalog(filename=catalog_filename)
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    #time_bins = np.linspace(0, 8, num=8*24*2 + 1) * u.day
    #time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-09-23 12:00')
    end_date = Time('2017-10-07 12:00')
    #date = start_date + time_bins

    date = time.compute_time(date_start=start_date, date_end=end_date, time_steps=15*u.min, location=location, only_night=True)
    objectives = np.ones(len(sources)) * date.shape[0] // len(sources)


    sun_intensity = sun.intensity(date=date, location=location)
    moon_intensity = moon.intensity(date=date, location=location)
    source_intensity = np.zeros((len(sources), date.shape[0]))
    altaz_moon = moon.compute_moon_position(date=date, location=location)

    for i, source in enumerate(sources):

        source_intensity[i] = gamma_source.intensity(date=date, location=location, ra=source['ra'], dec=source['dec'], altaz_moon=altaz_moon)

    sources_visibility = scheduler.compute_source_visibility(source_intensity, sun_intensity, moon_intensity)

    availability, schedule_quality = scheduler.find_quality_schedule(sources_visibility=sources_visibility)
    availability, schedule_quality_priority = scheduler.find_dynamic_priority_quality_schedule(
        sources_visibility=sources_visibility, objectives=objectives)

    writer.write_schedule(schedule=schedule_quality, sources=sources, dates=date, filename=output_directory + 'quality_schedule_%s.txt' % units_output, units=units_output)
    writer.write_schedule(schedule=schedule_quality_priority, sources=sources, dates=date, filename=output_directory + 'dynamic_priority_schedule_%s.txt' % units_output, units=units_output)

    fig_1, axis_1 = plt.subplots()
    axis_1.set_title('Quality schedule')
    fig_2, axis_2 = plt.subplots()
    axis_2.set_title('Priority schedule')
    fig_3, axis_3 = plt.subplots()
    axis_3.set_title('Dynamic quality-priority schedule')
    for i, source in enumerate(sources):

        if np.any(schedule_quality[i] > 0):
            axis_1.plot_date(date.plot_date, schedule_quality[i], label=source['name'], linestyle='None'
                             , marker='.')

        if np.any(schedule_quality_priority[i]>0):
            axis_3.plot_date(date.plot_date, schedule_quality_priority[i], label=source['name'], linestyle='None'
                             , marker='.')

    a = scheduler.compute_schedule_efficiency(sources_visibility=sources_visibility, schedule=schedule_quality)
    b = scheduler.compute_schedule_efficiency(sources_visibility=sources_visibility, schedule=schedule_quality_priority)

    print(a)
    print(b)

    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[occupancy]')
    axis_1.set_ylim([0.1, 1.2])
    axis_1.legend(loc='best')

    axis_2.set_xlabel('UTC time')
    axis_2.set_ylabel('[occupancy]')
    axis_2.set_ylim([0.1, 1.2])
    axis_2.legend(loc='best')

    axis_3.set_xlabel('UTC time')
    axis_3.set_ylabel('[occupancy]')
    axis_3.set_ylim([0.1, 1.2])
    axis_3.legend(loc='best')

    plt.gcf().autofmt_xdate()
    plt.show()
