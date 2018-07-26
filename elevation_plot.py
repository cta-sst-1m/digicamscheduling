import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon, sun, environement
from digicamscheduling.core.environement import interpolate_environmental_limits
from digicamscheduling.utils import time
from digicamscheduling.display.plot import plot_source_2d, plot_sun_2d, \
    plot_elevation
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.dates import date2num


def main(sources_filename, location_filename, environment_filename,
         start_date, end_date, time_steps):
    sources = reader.read_catalog(sources_filename)
    # sources = [sources[0]]
    coordinates = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates)

    alt_trees, az_trees = reader.read_environmental_limits(
        environment_filename)
    alt_trees = alt_trees * u.deg
    az_trees = az_trees * u.deg
    env_limits = interpolate_environmental_limits(alt_trees,
                                                               az_trees)

    start_date = Time(start_date)  # time should be 00:00
    end_date = Time(end_date)  # time should be 00:00
    hours = np.arange(0, 1, time_steps.to(u.day).value) * u.day
    hours = hours.to(u.hour)

    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_steps=time_steps, only_night=False)

    days = date.reshape(-1, len(hours))
    days = days.datetime
    days = date2num(days[:, 0])
    extent = [days.min(), days.max(), hours.value.min(), hours.value.max()]

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_elevation = moon_position.alt
    moon_phase = moon.compute_moon_phase(date=date)
    sun_position = sun.compute_sun_position(date=date, location=location)
    sun_elevation = sun_position.alt

    observability = (sun_elevation < -12 * u.deg) * np.cos(moon_elevation)
    observability *= (1 - moon_phase) * (moon_elevation < 0 * u.deg)

    for i, source in tqdm(enumerate(sources), total=len(sources),
                          desc='Source'):
        temp = gamma_source.compute_source_position(date=date,
                                                    location=location,
                                                    ra=source['ra'],
                                                    dec=source['dec'])
        source_elevation = temp.alt
        source_azimuth = temp.az
        is_above_trees = environement.is_above_environmental_limits(
            source_elevation, source_azimuth, env_limits)
        moon_separation = temp.separation(moon_position)

        source_visibility = is_above_trees * np.sin(source_elevation)
        source_visibility *= observability * (moon_separation > 10 * u.deg)

        source_elevation = source_elevation.reshape(-1, len(hours))
        moon_separation = moon_separation.reshape(-1, len(hours))
        source_visibility = source_visibility.reshape(-1, len(hours))
        # is_above_trees = is_above_trees.reshape(-1, len(hours))

        fig_1 = plt.figure()
        axes_1 = fig_1.add_subplot(111)
        fig_2 = plt.figure()
        axes_2 = fig_2.add_subplot(111)

        plot_source_2d(source_elevation, coordinates, source=source,
                       extent=extent, axes=axes_1,
                       vmin=-90, vmax=90, c_label='elevation [deg]')

        plot_source_2d(source_visibility, coordinates, source=source,
                       extent=extent, vmin=0, vmax=1, axes=axes_2,
                       c_label='visibility []')

        fig_1.savefig(source['name'] + '_elevation.png')
        fig_2.savefig(source['name'] + '_visibility.png')

        plt.show()

    observability = observability.reshape(-1, len(hours))
    moon_elevation = moon_elevation.reshape(-1, len(hours))
    moon_phase = moon_phase.reshape(-1, len(hours))
    sun_elevation = sun_elevation.reshape(-1, len(hours))

    plot_sun_2d(sun_elevation, coordinates, extent=extent)

    plot_source_2d(observability, coordinates, extent=extent,
                   c_label='Observability []', vmin=0, vmax=1)

    plot_source_2d(moon_elevation.value, coordinates, extent=extent,
                   vmin=-90, vmax=90, c_label='Moon elevation [deg]',
                   cmap=plt.get_cmap('RdYlGn_r'))
    plot_source_2d(moon_phase, coordinates, extent=extent,
                   c_label='Moon phase []', vmin=0, vmax=1,
                   cmap=plt.get_cmap('RdYlGn_r'))

    plt.show()


if __name__ == '__main__':
    start_date = '2018-01-01'
    end_date = '2018-12-31'
    time_step = 60 * u.minute

    location_filename = 'digicamscheduling/config/location_krakow.txt'
    sources_filename = 'digicamscheduling/config/catalog.txt'
    environment_filename = 'digicamscheduling/config/environmental_limitation.txt'

    main(location_filename=location_filename,
         sources_filename=sources_filename,
         environment_filename=environment_filename,
         start_date=start_date,
         end_date=end_date,
         time_steps=time_step)
