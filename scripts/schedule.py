import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.core.scheduler import compute_source_visibility, find_priority_schedule, compute_schedule_efficiency, find_quality_schedule, find_dynamic_priority_quality_schedule
import digicamscheduling.core.sun as sun
import digicamscheduling.core.moon as moon
import digicamscheduling.core.gamma_source as gamma_source
import matplotlib.pyplot as plt
from digicamscheduling.io import writer
from digicamscheduling.io import reader

if __name__ == '__main__':

    output_directory = '../digicamscheduling/samples/'

    catalog_filename = '../digicamscheduling/config/fact_catalog.txt'
    sources = reader.read_catalog(filename=catalog_filename)

    units_output = 'deg'

    weights = np.array([1, 3, 2])
    weights = weights / np.sum(weights)

    location_filename = '../digicamscheduling/config/location_krakow.txt'
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 8, num=8*24*2 + 1) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-08-26 12:00')

    sun_intensity = np.zeros(time_bins.shape)
    moon_intensity = np.zeros(time_bins.shape)
    source_intensity = np.zeros((len(sources), time_bins.shape[0]))

    for i, time in enumerate(time_bins):

        date = start_date + time
        sun_intensity[i] = sun.intensity(date=date, location=location)
        moon_intensity[i] = moon.intensity(date=date, location=location)

        for j, source in enumerate(sources):

            source_intensity[j, i] = gamma_source.intensity(date=date, location=location, ra=source['ra'], dec=source['dec'])

    sources_visibility = compute_source_visibility(source_intensity, sun_intensity, moon_intensity)
    objectives = np.ones(len(sources)) * 4 * u.hour

    print(objectives.to('day'))

    if np.sum(objectives) > (time_bins[-1] - time_bins[0]):

        print('Objectives are impossible')

    objectives = (objectives // time_interval)

    availability, schedule = find_priority_schedule(sources_visibility, objectives)

    fig, axis_1 = plt.subplots()
    # axis_1.plot_date((start_date + time_bins).plot_date, availability, label='schedule', linestyle='--', color='k', drawstyle='steps-mid', marker='None')

    for j, source in enumerate(sources):

        axis_1.plot_date((start_date + time_bins).plot_date, schedule[j], label=source['name'], linestyle='-', drawstyle='steps-mid', marker='None')

    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[occupancy]')
    axis_1.set_ylim([0, 1.2])
    plt.legend(loc='best')
    plt.gcf().autofmt_xdate()

    priority_schedule_efficiency = compute_schedule_efficiency(sources_visibility, schedule) * time_interval

    print(priority_schedule_efficiency)

    writer.write_schedule(schedule, sources, start_date, time_bins, output_directory + 'priority_schedule_%s.txt' % units_output, units=units_output)


    availability, schedule = find_quality_schedule(sources_visibility)

    fig, axis_1 = plt.subplots()
    # axis_1.plot_date((start_date + time_bins).plot_date, availability, label='schedule', linestyle='--', color='k', drawstyle='steps-mid', marker='None')

    for j, source in enumerate(sources):
        axis_1.plot_date((start_date + time_bins).plot_date, schedule[j], label=source['name'], linestyle='-',
                         drawstyle='steps-mid', marker='None')

    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[occupancy]')
    axis_1.set_ylim([0, 1.2])
    plt.legend(loc='best')
    plt.gcf().autofmt_xdate()

    quality_schedule_efficiency = compute_schedule_efficiency(sources_visibility, schedule) * time_interval

    print(quality_schedule_efficiency)

    writer.write_schedule(schedule, sources, start_date, time_bins, output_directory + 'quality_schedule_%s.txt' % units_output, units=units_output)

    availability, schedule = find_dynamic_priority_quality_schedule(sources_visibility, objectives)

    fig, axis_1 = plt.subplots()
    # axis_1.plot_date((start_date + time_bins).plot_date, availability, label='schedule', linestyle='--', color='k', drawstyle='steps-mid', marker='None')

    for j, source in enumerate(sources):
        axis_1.plot_date((start_date + time_bins).plot_date, schedule[j], label=source['name'], linestyle='-',
                         drawstyle='steps-mid', marker='None')

    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('[occupancy]')
    axis_1.set_ylim([0, 1.2])
    plt.legend(loc='best')
    plt.gcf().autofmt_xdate()

    dynamic_priority_schedule_efficiency = compute_schedule_efficiency(sources_visibility, schedule) * time_interval

    print(dynamic_priority_schedule_efficiency)

    writer.write_schedule(schedule, sources, start_date, time_bins, output_directory +'dynamic_priority_schedule_%s.txt' %units_output, units=units_output)


    plt.show()
