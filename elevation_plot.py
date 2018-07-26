import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon, sun, environement
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

    env_limits = environement.interpolate_environmental_limits(alt_trees,
                                                               az_trees)

    start_date = Time(start_date)  # time should be 00:00
    end_date = Time(end_date)  # time should be 00:00
    hours = np.arange(0, 1, time_steps.to(u.day).value) * u.day
    hours = hours.to(u.hour)

    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_steps=time_steps, only_night=False)

    source_elevation = np.zeros((len(sources), date.shape[0])) * u.deg
    source_azimuth = np.zeros((len(sources), date.shape[0])) * u.deg
    is_above_trees = np.zeros((len(sources), date.shape[0]), dtype=bool)

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_elevation = moon_position.alt
    moon_phase = moon.compute_moon_phase(date=date)
    sun_position = sun.compute_sun_position(date=date, location=location)
    sun_elevation = sun_position.alt

    moon_separation = np.zeros((len(sources), len(date))) * u.deg

    for i, source in tqdm(enumerate(sources), total=len(sources),
                          desc='Source'):
        temp = gamma_source.compute_source_position(date=date,
                                                    location=location,
                                                    ra=source['ra'],
                                                    dec=source['dec'])
        source_elevation[i] = temp.alt
        source_azimuth[i] = temp.az
        is_above_trees[i] = environement.is_above_environmental_limits(
            temp.alt, temp.az, env_limits)
        moon_separation[i] = temp.separation(moon_position)

    source_elevation = source_elevation.reshape(len(sources), -1, len(hours))
    moon_separation = moon_separation.reshape(len(sources), -1, len(hours))
    is_above_trees = is_above_trees.reshape(len(sources), -1, len(hours))
    moon_elevation = moon_elevation.reshape(-1, len(hours))
    moon_phase = moon_phase.reshape(-1, len(hours))
    sun_elevation = sun_elevation.reshape(-1, len(hours))

    observability = (sun_elevation < -12 * u.deg) * np.cos(moon_elevation) \
                    * (1 - moon_phase) * (moon_elevation < 0 * u.deg)

    source_visibility = is_above_trees * np.sin(
        source_elevation) * observability

    source_visibility = source_visibility * (moon_separation > 10 * u.deg)

    date = date.reshape(-1, len(hours))
    date = date.datetime
    days = date2num(date[:, 0])

    extent = [days.min(), days.max(),
              hours.value.min(), hours.value.max()]

    plot_sun_2d(sun_elevation, coordinates, extent=extent)

    plot_source_2d(observability, coordinates, extent=extent,
                   c_label='Observability []', vmin=0, vmax=1)

    plot_source_2d(moon_elevation.value, coordinates, extent=extent,
                   vmin=-90, vmax=90, c_label='Moon elevation [deg]',
                   cmap=plt.get_cmap('RdYlGn_r'))
    plot_source_2d(moon_phase, coordinates, extent=extent,
                   c_label='Moon phase []', vmin=0, vmax=1,
                   cmap=plt.get_cmap('RdYlGn_r'))

    for i, source in enumerate(sources):
        # az = source_azimuth[i]
        alt = source_elevation[i]
        moon_sep = moon_separation[i]
        visibility = source_visibility[i]

        plot_source_2d(alt, coordinates, source=source, extent=extent,
                       vmin=-90, vmax=90, c_label='elevation [deg]')

        plot_source_2d(visibility, coordinates, source=source,
                       extent=extent, vmin=0, vmax=1,
                       c_label='visibility []')

        plot_source_2d(moon_sep, coordinates, source=source, extent=extent,
                       vmin=0, vmax=180, c_label='Moon separation [deg]')

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
