import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from core.scheduler import compute_source_visibility, find_priority_schedule, compute_schedule_efficiency, find_quality_schedule, find_dynamic_priority_quality_schedule
import core.sun as sun
import core.moon as moon
import core.gamma_source as gamma_source
import matplotlib.pyplot as plt
from output.writer import write_schedule


if __name__ == '__main__':

    sources = [{'ra': 11.074266 * u.deg, 'dec': 38.208801 * u.deg, 'name': 'Mrk 421'},
               {'ra': 5.575539 * u.deg, 'dec': 22.014500 * u.deg, 'name': 'Crab'},
               {'ra': 16.897867 * u.deg, 'dec': 39.760201 * u.deg, 'name': 'Mrk 501'}]

    weights = np.array([1, 3, 2])
    weights = weights / np.sum(weights)

    coordinates_krakow = {'lat': 50.090815 * u.deg, 'lon': 19.887937 * u.deg, 'height': 214.034 * u.m}
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 7, num=7*24*2 + 1) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-08-27 00:00')

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
    objectives = np.array([10, 10, 20]) * u.hour

    print(objectives)

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

    write_schedule(schedule, sources, start_date, time_bins, 'priority_schedule.txt')


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

    write_schedule(schedule, sources, start_date, time_bins, 'quality_schedule.txt')

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

    write_schedule(schedule, sources, start_date, time_bins, 'dynamic_priority_schedule.txt')


    plt.show()
