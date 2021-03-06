"""
Plot catalog

Usage:
  digicamscheduling-catalog [options]

Options:
    --help                       Show this screen.
    --start_date=DATE            Starting date (UTC) YYYY-MM-DD
                                 (Should start at midnight !)
                                 [default: 2018-01-01 00:00]
    --end_date=DATE              Ending date (UTC) YYYY-MM-DD
                                 (Should end at midnight !)
                                 [default: 2019-01-01 00:00]
    --time_step=MINUTES          Time steps in minutes
                                 [default: 60]
    --output_path=PATH           Path to save the figure
    --location_filename=PATH     PATH for location config file
    --sources_filename=PATH      PATH for catalog
    --environment_filename=PATH  PATH for environmental limitations
    --hide                       Hide the plots
"""
from docopt import docopt
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon, sun, environement
from digicamscheduling.core.environement import \
    interpolate_environmental_limits, compute_observability
from digicamscheduling.utils import time
from digicamscheduling.utils.docopt import convert_commandline_arguments
from digicamscheduling.display.plot import plot_source_2d
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.dates import date2num
import os


def main(sources_filename, location_filename, environment_filename,
         start_date, end_date, time_step, output_path, hide=False):

    sources = reader.read_catalog(sources_filename)
    coordinates = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates)

    alt_trees, az_trees = reader.read_environmental_limits(
        environment_filename)
    alt_trees = alt_trees * u.deg
    az_trees = az_trees * u.deg
    env_limits = interpolate_environmental_limits(alt_trees, az_trees)
    start_date = Time(start_date)  # time should be 00:00
    end_date = Time(end_date)  # time should be 00:00
    hours = np.arange(0, 1, time_step.to(u.day).value) * u.day
    hours = hours.to(u.hour)

    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_step=time_step, only_night=False)

    days = date.reshape(-1, len(hours))
    days = days.datetime
    days = date2num(days[:, 0])
    extent = [days.min(), days.max(), hours.value.min(), hours.value.max()]

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_elevation = moon_position.alt
    moon_phase = moon.compute_moon_phase(date=date)
    sun_position = sun.compute_sun_position(date=date, location=location)
    sun_elevation = sun_position.alt

    observability = compute_observability(sun_elevation, moon_elevation,
                                          moon_phase, use_moon=True)

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
        source_visibility = source_visibility.reshape(-1, len(hours))

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

        if output_path is not None:

            filename = os.path.join(output_path, source['name'])

            fig_1.savefig(filename + '_elevation.png')
            fig_2.savefig(filename + '_visibility.png')

    if not hide:

        plt.show()


def entry():

    kwargs = docopt(__doc__)
    kwargs = convert_commandline_arguments(kwargs)
    main(**kwargs)
