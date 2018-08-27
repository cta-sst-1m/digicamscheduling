"""
Plot catalog

Usage:
  digicamscheduling-elevation [options]

Options:
 -h --help                   Show this screen.
 --start_date=DATE            Starting date (UTC) YYYY-MM-DD HH:MM:SS
                              [default: 2018-01-01 00:00:00]
 --end_date=DATE              Ending date (UTC) YYYY-MM-DD HH:MM:SS
                              [default: 2018-01-15 00:00:00]
 --time_step=MINUTES          Time steps in minutes
                              [default: 5]
 --output_path=PATH           Path to save the figure. If not specified the
                              figures will not be saved
 --location_filename=PATH     PATH for location config file
                              [default: digicamscheduling/config/location_krakow.txt]
 --sources_filename=PATH      PATH for catalog
                              [default: digicamscheduling/config/catalog.json]
 --environment_filename=PATH  PATH for environmental limitations
                              [default: digicamscheduling/config/environmental_limitation.txt]
 --show                       View directly the plot
 --threshold=N                Threshold for visibility
                              [default: 0.0]
 --use_moon                   Choose to use moon to compute source visibility
                              [Default: False]
"""
from docopt import docopt
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon, sun
from digicamscheduling.core.environement import interpolate_environmental_limits, is_above_environmental_limits, compute_observability
from digicamscheduling.utils import time
from digicamscheduling.display.plot import plot_elevation, plot_source
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from tqdm import tqdm
import os


def main(sources_filename, location_filename, environment_filename,
         start_date, end_date, time_steps, output_path, use_moon, show=False,
         threshold=0.5):

    sources = reader.read_catalog(sources_filename)
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

    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_steps=time_steps, location=location,
                             only_night=True)

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_elevation = moon_position.alt
    moon_phase = moon.compute_moon_phase(date=date)
    sun_position = sun.compute_sun_position(date=date, location=location)
    sun_elevation = sun_position.alt

    observability = compute_observability(sun_elevation, moon_elevation,
                                          moon_phase, use_moon=use_moon)

    fig_1 = plt.figure()
    axes_1 = fig_1.add_subplot(111)
    fig_2 = plt.figure()
    axes_2 = fig_2.add_subplot(111)

    color = iter(cm.rainbow(np.linspace(0, 1, num=len(sources))))

    for i, source in tqdm(enumerate(sources), total=len(sources),
                          desc='Source'):

        c = next(color)

        temp = gamma_source.compute_source_position(date=date,
                                                    location=location,
                                                    ra=source['ra'],
                                                    dec=source['dec'])
        source_elevation = temp.alt
        source_azimuth = temp.az
        is_above_trees = is_above_environmental_limits(
            source_elevation, source_azimuth, env_limits)
        moon_separation = temp.separation(moon_position)

        source_visibility = is_above_trees * np.sin(source_elevation)
        source_visibility *= observability * (moon_separation > 10 * u.deg)
        source_visibility *= source['weight']

        if np.any(source_visibility >= threshold):

            label = source['name']
            plot_elevation(date, source_elevation, axes=axes_1, color=c,
                           label=label)
            plot_source(date, source_visibility,
                        axes=axes_2, color=c, y_label='visibility []',
                        ylim=[threshold, 1], label=label)

    if output_path is not None:

        fig_1.savefig(os.path.join(output_path, 'elevation.png'))
        fig_2.savefig(os.path.join(output_path, 'visibility.png'))

    if show:

        plt.show()


def entry():

    args = docopt(__doc__)

    main(location_filename=args['--location_filename'],
         sources_filename=args['--sources_filename'],
         environment_filename=args['--environment_filename'],
         start_date=args['--start_date'],
         end_date=args['--end_date'],
         time_steps=float(args['--time_step']) * u.minute,
         output_path=args['--output_path'],
         show=args['--show'],
         threshold=float(args['--threshold']),
         use_moon=args['--use_moon'])


if __name__ == '__main__':

    start_date = '2018-06-26'
    end_date = '2018-07-10'
    time_step = 1 * u.minute
    output_path = 'figures/'
    show = False

    location_filename = 'digicamscheduling/config/location_krakow.txt'
    sources_filename = 'digicamscheduling/config/catalog.json'
    environment_filename = 'digicamscheduling/config/' \
                           'environmental_limitation.txt'

    main(location_filename=location_filename,
         sources_filename=sources_filename,
         environment_filename=environment_filename,
         start_date=start_date,
         end_date=end_date,
         time_steps=time_step,
         output_path=output_path,
         show=show)
