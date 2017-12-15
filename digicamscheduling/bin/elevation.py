#!/usr/bin/env python

import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.io import reader
from digicamscheduling.core import gamma_source, moon
from digicamscheduling.utils import time
import digicamscheduling.display.plot as display
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


def plot_elevation():

    sources_filename = 'digicamscheduling/config/' + 'catalog.txt'
    location_filename = 'digicamscheduling/config/' + 'location_krakow.txt'

    sources = reader.read_catalog(sources_filename)

    sources_names = [source['name'] for source in sources]
    sources = [sources[sources_names.index(
        'Crab')]]  # , sources[sources_names.index('3FGL J0509.4+0541')]]# , sources[sources_names.index('1ES 1959+650')]]

    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    start_date = Time('2017-11-27 12:00')
    end_date = Time('2017-11-28 12:00')
    time_steps = 1 * u.min
    date = time.compute_time(date_start=start_date, date_end=end_date, time_steps=time_steps, location=location)
    source_elevation = np.zeros((len(sources), date.shape[0])) * u.deg
    source_azimut = np.zeros((len(sources), date.shape[0])) * u.deg

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_phase = moon.compute_moon_phase(date=date)
    moon_phase = np.mean(moon_phase)

    for i, source in enumerate(sources):
        temp = gamma_source.compute_source_position(date=date, location=location, ra=source['ra'], dec=source['dec'])
        source_elevation[i] = temp.alt
        source_azimut[i] = temp.az

    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()

    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')

    color = iter(cm.rainbow(np.linspace(0, 1, len(sources) + 1)))
    for i, source in enumerate(sources):

        c = next(color)

        if np.any(source_elevation[i] > 45 * u.deg) or source['name'] == 'Crab':
            alt = source_elevation[i]
            az = source_azimut[i]
            display.plot_azimuth(date, az, axis=axis_1,
                                 label=source['name'] + ', flux : {} [c.u.]'.format(source['flux']), color=c)
            display.plot_elevation(date, alt, axis=axis_2,
                                   label=source['name'] + ', flux : {} [c.u.]'.format(source['flux']), color=c)
            display.plot_trajectory(az, alt, axis=axis_3,
                                    label=source['name'] + ', flux : {} [c.u.]'.format(source['flux']), color=c)

        if i == 0:
            c = next(color)
            display.plot_elevation(date, moon_position.alt, axis=axis_2,
                                   label='Moon, phase : {:0.1f}'.format(moon_phase), color='k')
            display.plot_azimuth(date, moon_position.az, axis=axis_1, label='Moon, phase : {:0.1f}'.format(moon_phase),
                                 color='k')
            display.plot_trajectory(moon_position.az, moon_position.alt, axis=axis_3,
                                    label='Moon, phase : {:0.1f}'.format(moon_phase), color='k')

    plt.show()


if __name__ == '__main__':

    plot_elevation()
