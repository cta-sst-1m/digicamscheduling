"""
Plot catalog

Usage:
  digicamscheduling-elevation [options]

Options:
 -h --help                   Show this screen.
 --start_date=DATE            Starting date (UTC) YYYY-MM-DD HH:MM:SS
                              [default: 2018-01-01 00:00:00]
 --end_date=DATE              Ending date (UTC) YYYY-MM-DD HH:MM:SS
                              [default: 2018-12-31 00:00:00]
 --time_step=MINUTES          Time steps in minutes
                              [default: 60]
 --output_path=PATH           Path to save the figure
 --location_filename=PATH     PATH for location config file
 --hide                       Hide the plot
"""
from docopt import docopt
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import moon, sun
from digicamscheduling.core.environement import compute_observability
from digicamscheduling.utils import time
from digicamscheduling.display.plot import plot_source_2d, plot_sun_2d
from digicamscheduling.utils.docopt import convert_commandline_arguments
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import os


def main(location_filename, start_date, end_date, time_step, output_path,
         hide=False):

    coordinates = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates)

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
                                          moon_phase)

    observability = observability.reshape(-1, len(hours))
    moon_elevation = moon_elevation.reshape(-1, len(hours))
    moon_phase = moon_phase.reshape(-1, len(hours))
    sun_elevation = sun_elevation.reshape(-1, len(hours))

    fig_1 = plt.figure()
    axes_1 = fig_1.add_subplot(111)
    fig_2 = plt.figure()
    axes_2 = fig_2.add_subplot(111)
    fig_3 = plt.figure()
    axes_3 = fig_3.add_subplot(111)
    fig_4 = plt.figure()
    axes_4 = fig_4.add_subplot(111)

    plot_sun_2d(sun_elevation, coordinates, extent=extent, axes=axes_1)

    plot_source_2d(observability, coordinates, extent=extent,
                   c_label='Observability []', vmin=0, vmax=1, axes=axes_2)

    plot_source_2d(moon_elevation.value, coordinates, extent=extent,
                   vmin=-90, vmax=90, c_label='Moon elevation [deg]',
                   cmap=plt.get_cmap('RdYlGn_r'), axes=axes_3)
    plot_source_2d(moon_phase, coordinates, extent=extent,
                   c_label='Moon phase []', vmin=0, vmax=1,
                   cmap=plt.get_cmap('RdYlGn_r'), axes=axes_4)

    if output_path is not None:

        fig_1.savefig(os.path.join(output_path, 'sun_elevation.png'))
        fig_2.savefig(os.path.join(output_path, 'observability.png'))
        fig_3.savefig(os.path.join(output_path, 'moon_elevation.png'))
        fig_4.savefig(os.path.join(output_path, 'moon_phase.png'))

    if not hide:

        plt.show()


def entry():

    kwargs = docopt(__doc__)
    kwargs = convert_commandline_arguments(kwargs)
    main(**kwargs)
