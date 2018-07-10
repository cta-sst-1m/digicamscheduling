import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon, sun, environement
from digicamscheduling.utils import time
import digicamscheduling.display.plot as display
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from tqdm import tqdm
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange, date2num


def plot_source(source_elevation, source_name, coordinates, extent=None):

    cmap = plt.get_cmap('RdYlGn')
    cmap.set_bad(color='k', alpha=1.)
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.set_title('Site : ({}, {}, {}) [lat, lon, alt]'
                   '\nSource : {}'.format(coordinates['lat'],
                                          coordinates['lon'],
                                          coordinates['height'],
                                          source_name), fontsize=12)

    ax = axes.imshow(source_elevation.T, aspect='auto', origin='lower',
                     extent=extent, vmin=-90, vmax=90, cmap=cmap)

    axes.xaxis_date()
    axes.xaxis.set_minor_formatter(DateFormatter('%Y-%m-%d'))
    axes.set_xlabel('date')
    axes.set_ylabel('hour [UTC]')
    fig.autofmt_xdate()
    fig.colorbar(ax, label='elevation [deg]', extend='both')

    return fig, axes


def main(sources_filename='digicamscheduling/config/catalog.txt',
         location_filename='digicamscheduling/config/location_krakow.txt',
         environment_filename='digicamscheduling/config/environmental_limitation.txt'):

    sources = reader.read_catalog(sources_filename)
    # sources = [sources[4]]
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
    time_steps = 30 * u.minute
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

    mask_sun_and_moon = (sun_elevation < -18 * u.deg) *  \
                                        (moon_elevation < 0. * u.deg)

    date = date.reshape(-1, len(hours))
    date = date.datetime
    days = date2num(date[:, 0])

    extent = [days.min(), days.max(),
              hours.value.min(), hours.value.max()]

    plot_source(sun_elevation.value, 'Sun', coordinates, extent=extent)
    plot_source(moon_elevation.value, 'Moon', coordinates, extent=extent)

    for i, source in enumerate(sources):

        az = source_azimuth[i].value
        alt = source_elevation[i].value
        mask = is_above_trees[i]
        mask = mask * mask_sun_and_moon

        alt = np.ma.masked_array(alt, mask=~mask)

        plot_source(alt, source['name'],
                    coordinates, extent=extent)

        display.plot_trajectory(source_elevation[i].ravel(), source_azimuth[i].ravel())

    plt.show()


if __name__ == '__main__':
    main(location_filename='digicamscheduling/config/location_krakow.txt')
    # main()
