import numpy as np
import os
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon, sun, environement
from digicamscheduling.utils import time
import digicamscheduling.display.plot as display
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.dates import DateFormatter, date2num


def plot_source(source_elevation, coordinates, source=None,
                c_label='elevation [deg]', extent=None, **kwargs):

    cmap = plt.get_cmap('RdYlGn')
    cmap.set_bad(color='k', alpha=1.)
    fig = plt.figure()
    axes = fig.add_subplot(111)

    title = 'Site : ({}, {}, {}) [lat, lon, alt] deg'.format(
        coordinates['lat'], coordinates['lon'], coordinates['height']
    )

    if source is not None:

        title += '\nSource : {} ({}, {}) [ra, dec] deg'.format(source['name'],
                                                               source['ra'],
                                                                source['dec'])

    axes.set_title(title, fontsize=20)

    ax = axes.imshow(source_elevation.T, aspect='auto', origin='lower',
                     extent=extent, cmap=cmap, **kwargs)

    axes.xaxis_date()
    axes.xaxis.set_minor_formatter(DateFormatter('%Y-%m-%d'))
    axes.set_xlabel('date')
    axes.set_ylabel('hour [UTC]')
    axes.set_yticks(np.arange(0, 24, 1))
    fig.autofmt_xdate()
    fig.colorbar(ax, label=c_label, extend='both')

    return fig, axes


def plot_sun(sun_elevation, coordinates, extent, **kwargs):

    fig, axes = plot_source(sun_elevation, coordinates, extent=extent,
                            vmin=-90, vmax=90, **kwargs,
                            c_label='Sun elevation [deg]')

    cs = axes.contour(sun_elevation.T, levels=[-18., -12., -6., -0.],
                      extent=extent, cmap='binary_r')
    sun_max = np.argmax(sun_elevation, axis=1)
    days = np.linspace(extent[0], extent[1], num=len(sun_elevation))
    hours = np.linspace(extent[2], extent[3], num=sun_elevation.shape[1])
    hours = hours[sun_max]

    axes.plot(days, hours, color='r')

    contour_labels = ['Astronomical', 'Nautical', 'Civil',
                      'Horizon']
    fmt = dict(zip(cs.levels, contour_labels))
    axes.clabel(cs, fmt=fmt)

    return fig, axes


def main(sources_filename='digicamscheduling/config/catalog.txt',
         location_filename='digicamscheduling/config/location_krakow.txt',
         environment_filename='digicamscheduling/config/environmental_limitation.txt'):

    sources = reader.read_catalog(sources_filename)
    # sources = [sources[-1]]
    coordinates = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates)

    alt_trees, az_trees = reader.read_environmental_limits(
        environment_filename)
    alt_trees = alt_trees * u.deg
    az_trees = az_trees * u.deg

    env_limits = environement.interpolate_environmental_limits(alt_trees,
                                                               az_trees)

    start_date = Time('2018-01-01') # time should be 00:00
    end_date = Time('2018-12-31')  # time should be 00:00
    time_steps = 15 * u.minute
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

    source_elevation = source_elevation.reshape(len(sources), -1, len(hours))
    source_azimuth = source_azimuth.reshape(len(sources), -1, len(hours))
    is_above_trees = is_above_trees.reshape(len(sources), -1, len(hours))
    moon_elevation = moon_elevation.reshape(-1, len(hours))
    moon_phase = moon_phase.reshape(-1, len(hours))
    sun_elevation = sun_elevation.reshape(-1, len(hours))

    observability = (sun_elevation < -12 * u.deg) * np.cos(moon_elevation) \
                    * (1 - moon_phase)

    source_visibility = is_above_trees * np.sin(source_elevation) * observability

    date = date.reshape(-1, len(hours))
    date = date.datetime
    days = date2num(date[:, 0])

    extent = [days.min(), days.max(),
              hours.value.min(), hours.value.max()]

    sun_fig, sun_axes = plot_sun(sun_elevation.value, coordinates,
                                 extent=extent)

    plot_source(observability, coordinates, extent=extent,
                c_label='Observability []', vmin=0, vmax=1)

    plot_source(moon_elevation.value, coordinates, extent=extent,
                vmin=-90, vmax=90)
    plot_source(moon_phase, coordinates, extent=extent,
                c_label='Moon phase []', vmin=0, vmax=1)

    for i, source in enumerate(sources):

        az = source_azimuth[i].value
        alt = source_elevation[i].value
        visibility = source_visibility[i]

        fig_1, temp = plot_source(alt, coordinates, source=source,
                                  extent=extent, vmin=-90, vmax=90)

        fig_2, temp = plot_source(visibility, coordinates, source=source,
                                  extent=extent, vmin=0, vmax=1,
                                  c_label='visibility []')

        # path = '/home/alispach/figures/visibility/2018_observations/sources'
        # file_1 = os.path.join(path, source['name'] + '_elevation.svg')
        # file_2 = os.path.join(path, source['name'] + '_visibility.svg')
        # fig_1.savefig(file_1)
        # fig_2.savefig(file_2)

        # display.plot_trajectory(source_elevation[i].ravel(), source_azimuth[i].ravel())

    plt.show()


if __name__ == '__main__':
    main(location_filename='digicamscheduling/config/location_krakow.txt')
    # main()
