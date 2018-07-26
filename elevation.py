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
from tqdm import tqdm


def main(sources_filename='digicamscheduling/config/catalog.txt',
         location_filename='digicamscheduling/config/location_krakow.txt'):

    sources = reader.read_catalog(sources_filename)
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    start_date = Time('2018-07-27 00:00') # time should be 00:00
    end_date = Time('2018-08-03 00:00') # time should be 00:00
    time_steps = 1 * u.minute
    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_steps=time_steps, location=location)

    source_elevation = np.zeros((len(sources), date.shape[0])) * u.deg
    source_azimuth = np.zeros((len(sources), date.shape[0])) * u.deg

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_phase = moon.compute_moon_phase(date=date)

    for i, source in tqdm(enumerate(sources), total=len(sources),
                          desc='Source'):

        temp = gamma_source.compute_source_position(date=date,
                                                    location=location,
                                                    ra=source['ra'],
                                                    dec=source['dec'])
        source_elevation[i] = temp.alt
        source_azimuth[i] = temp.az

    color = iter(cm.rainbow(np.linspace(0, 1, len(sources))))
    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()

    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')

    for i, source in enumerate(sources):

        c = next(color)
        label = source['name'] + ', flux : {} [c.u.]'.format(source['flux'])

        if np.any(source_elevation[i] > 45 * u.deg) or source['name'] == 'Crab':
            alt = source_elevation[i]
            az = source_azimuth[i]

            display.plot_azimuth(date, az, axis=axis_1,
                                 label=label, color=c)
            display.plot_elevation(date, alt, axis=axis_2,
                                   label=label, color=c, linestyle='-',
                                   marker='None')
            display.plot_trajectory(az, alt, axes=axis_3,
                                    label=label, color=c)

        if i == 0:

            label_moon = 'Moon'
            display.plot_elevation(date, moon_position.alt, axis=axis_2,
                                   label=label_moon, color='k')
            display.plot_azimuth(date, moon_position.az, axis=axis_1,
                                 label=label_moon, color='k')
            display.plot_trajectory(moon_position.az, moon_position.alt,
                                    axes=axis_3, label=label_moon, color='k')

    plt.show()


if __name__ == '__main__':

    main(location_filename='digicamscheduling/config/location_krakow.txt')
    # main()
