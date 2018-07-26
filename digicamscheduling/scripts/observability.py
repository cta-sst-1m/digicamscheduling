"""
Plot catalog

Usage:
  digicamscheduling-elevation [options]

Options:
 -h --help                   Show this screen.
 --start_date=DATE            Starting date YYYY-MM-DD-HH:MM:SS
                              [default: 2018-01-01 00:00:00]
 --end_date=DATE              Ending date YYYY-MM-DD-HH:MM:SS
                              [default: 2018-12-31 00:00:00]
 --time_step=MINUTES          Time steps in minutes
                              [default: 60]
 --output_path=PATH           Path to save the figure
                              [default: .]
 --location_filename=PATH     PATH for location config file
                              [default: digicamscheduling/config/location_krakow.txt]
 --sources_filename=PATH      PATH for catalog
 --show                       View directly the plot
                              [default: False]
"""
from docopt import docopt
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import moon, sun
from digicamscheduling.utils import time
from digicamscheduling.display.plot import plot_source_2d, plot_sun_2d, plot_elevation
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import os


def main(location_filename, start_date, end_date, time_steps, output_path,
         show=False):

    coordinates = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates)

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

    fig_1.savefig(os.path.join(output_path, 'sun_elevation.png'))
    fig_2.savefig(os.path.join(output_path, 'observability.png'))
    fig_3.savefig(os.path.join(output_path, 'moon_elevation.png'))
    fig_4.savefig(os.path.join(output_path, 'moon_phase.png'))

    if show:

        plt.show()


def entry():

    args = docopt(__doc__)

    main(location_filename=args['--location_filename'],
         start_date=args['--start_date'],
         end_date=args['--end_date'],
         time_steps=float(args['--time_step']) * u.minute,
         output_path=args['--output_path'],
         show=args['--show'])


if __name__ == '__main__':
    start_date = '2018-08-01'
    end_date = '2018-08-31'
    time_step = 15 * u.minute
    output_path = 'figures/'
    show = True

    location_filename = 'digicamscheduling/config/location_krakow.txt'

    main(location_filename=location_filename,
         start_date=start_date,
         end_date=end_date,
         time_steps=time_step,
         output_path=output_path,
         show=show,)
