import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time


def compute_source_visibility(source_intensity, sun_intensity, moon_intensity, flux=None):

    source_visibility = source_intensity * (1. - sun_intensity) * (1. - moon_intensity)

    if flux is None:

        return source_visibility

    else:

        return source_visibility * flux


def compute_schedule_efficiency(sources_visibility, schedule, weights=None):

    if weights is None:

        observation_efficiency = np.sum(sources_visibility*schedule, axis=1)

    return observation_efficiency


def find_priority_schedule(sources_visibility, objectives):  # objectives should be in number of time intervalls

    availability = np.zeros(sources_visibility.shape[-1])
    schedule = np.zeros(sources_visibility.shape)
    what_to_observe = np.argsort(objectives)

    for source_id in reversed(what_to_observe):  # Loop over the sources in priority order

        n_observations = objectives[source_id]
        source_visibility = sources_visibility[source_id]

        if np.max(source_visibility) < 0.2:

            break

        when_to_observe = np.argsort(source_visibility)

        for time_bin in reversed(when_to_observe):

            if availability[time_bin] == 0:
                schedule[source_id][time_bin] = 1
                availability[time_bin] = 1

            if np.sum(schedule[source_id]) >= n_observations:

                break

    return availability.astype(bool), schedule.astype(bool)


def find_dynamic_priority_quality_schedule(sources_visibility, objectives, slew_time=0.2):  # objectives should be in number of time intervalls

    availability = np.zeros(sources_visibility.shape[-1])
    schedule = np.zeros(sources_visibility.shape)
    achievement = np.zeros(schedule.shape[0])

    current_source = 0
    slew_times = np.ones(schedule.shape[0]) * slew_time
    slew_times[current_source] = 0

    for time_bin in range(availability.shape[0]):

        if np.any(sources_visibility[..., time_bin] > 0.2):

            remaining_time = (objectives - achievement) / objectives
            priority_source = np.argmax(remaining_time * sources_visibility[..., time_bin] - slew_times)
            slew_times[current_source] = slew_time
            current_source = priority_source
            slew_times[current_source] = 0
            achievement[priority_source] += sources_visibility[priority_source, time_bin]
            schedule[priority_source, time_bin] = 1
            availability[time_bin] = 1

    return availability.astype(bool), schedule.astype(bool)


def find_quality_schedule(sources_visibility):

    availability = np.zeros(sources_visibility.shape[-1])
    schedule = np.zeros(sources_visibility.shape)
    previous_source = None

    for time_bin in range(sources_visibility.shape[1]):

        current_source = np.argmax(sources_visibility[:, time_bin])

        if current_source == previous_source:

            temp = np.argsort(sources_visibility[:, time_bin])
            new_source = temp[-2]
            if sources_visibility[new_source, time_bin] > 0.6:

                current_source = new_source

        if sources_visibility[current_source, time_bin] > 0.2:

            schedule[current_source, time_bin] = 1
            availability[time_bin] = 1

        previous_source = current_source

    return availability.astype(bool), schedule.astype(bool)


if __name__ == '__main__':

    time_bins = np.linspace(0, 7, num=10) * u.day
    sources = [{'ra': 11.074266 * u.deg, 'dec': 38.208801 * u.deg, 'name': 'Mrk 421'},
                          {'ra': 5.575539 * u.deg, 'dec': 22.014500 * u.deg, 'name': 'Crab'},
                          {'ra': 16.897867 * u.deg, 'dec':  39.760201 * u.deg, 'name': 'Mrk 501'}]

    weights = np.array([1, 3, 2])
    weights = weights / np.sum(weights)

    coordinates_krakow = {'lat': 50.090763 * u.deg, 'lon': 19.887956 * u.deg, 'height': 230 * u.m}
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 2, num=20) * u.day
    start_date = Time('2017-08-27 00:00')
    objectives = np.array([2., 6., 1.]) * u.hour


    # availability, best_schedule = find_schedule(location, start_date, time_bins, sources, weights, objectives, type='priority')


    # print(availability)
    # print(best_schedule)