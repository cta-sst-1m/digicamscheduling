"""
Plot Moon

Usage:
  digicamscheduling-moon [options]

Options:
 -h --help                   Show this screen.
 --start_date=DATE            Starting date (UTC) YYYY-MM-DD HH:MM:SS
 --end_date=DATE              Ending date (UTC) YYYY-MM-DD HH:MM:SS
 --time_step=MINUTES          Time steps in minutes
                              [default: 5]
 --output_path=PATH           Path to save the figure. If not specified the
                              figures will not be saved
 --location_filename=PATH     PATH for location config file
                              [default: digicamscheduling/config/location_krakow.txt]
 --show=BOOL                  View directly the plot
                              [default: 1]
"""
import os
from docopt import docopt
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.core import moon
from digicamscheduling.io import reader
import digicamscheduling.display.plot as display
import matplotlib.pyplot as plt

from digicamscheduling.utils import time
from digicamscheduling.utils.docopt import convert_commandline_arguments

def main(location_filename, start_date, end_date, time_step, output_path,
         show=False):

    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)
    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_step=time_step, location=location,
                             only_night=True)
    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_phase = moon.compute_moon_phase(date=date)

    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()
    fig_4 = plt.figure()
    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')
    axis_4 = fig_4.add_subplot(111)

    display.plot_azimuth(date, moon_position.az, axes=axis_1,
                         label='Moon')
    display.plot_elevation(date, moon_position.alt, axes=axis_2,
                           label='Moon')
    display.plot_trajectory(moon_position.az, moon_position.alt, axes=axis_3,
                            label='Moon')
    display.plot_moon_phase(date=date, phase=moon_phase, axes=axis_4,
                            label='Moon')

    if show:

        plt.show()

    if output_path is not None:

        fig_1.savefig(os.path.join(output_path, 'moon_azimuth.png'))
        fig_2.savefig(os.path.join(output_path, 'moon_elevation.png'))
        fig_3.savefig(os.path.join(output_path, 'moon_trajectory.png'))
        fig_4.savefig(os.path.join(output_path, 'moon_phase.png'))


def entry():

    kwargs = docopt(__doc__)
    print(kwargs)

    kwargs = convert_commandline_arguments(kwargs)

    main(**kwargs)


if __name__ == '__main__':

    location_filename = 'digicamscheduling/config/' + 'location_krakow.txt'
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    start_date = Time('2018-07-27 12:00')
    end_date = Time('2018-07-28 12:00')
    time_steps = 1 * u.min
    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_step=time_steps, location=location,
                             only_night=True)
    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_phase = moon.compute_moon_phase(date=date)

    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()
    fig_4 = plt.figure()
    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')
    axis_4 = fig_4.add_subplot(111)

    display.plot_azimuth(date, moon_position.az, axes=axis_1,
                         label='Moon')
    display.plot_elevation(date, moon_position.alt, axes=axis_2,
                           label='Moon')
    display.plot_trajectory(moon_position.az, moon_position.alt, axes=axis_3,
                            label='Moon')
    display.plot_moon_phase(date=date, phase=moon_phase, axes=axis_4,
                            label='Moon')

    plt.show()
