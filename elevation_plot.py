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
    sources = [sources[0]]
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    start_date = Time('2018-07-09 00:00') # time should be 00:00
    end_date = Time('2019-01-09 00:00')  # time should be 00:00
    time_steps = 1 * u.hour
    hours = np.arange(0, 24, 1) * u.hour

    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_steps=time_steps, only_night=False)

    print(date.shape)

    source_elevation = np.zeros((len(sources), date.shape[0])) * u.deg

    moon_position = moon.compute_moon_position(date=date, location=location)
    moon_phase = moon.compute_moon_phase(date=date)

    for i, source in tqdm(enumerate(sources), total=len(sources),
                          desc='Source'):
        temp = gamma_source.compute_source_position(date=date,
                                                    location=location,
                                                    ra=source['ra'],
                                                    dec=source['dec'])
        source_elevation[i] = temp.alt

    source_elevation = source_elevation.reshape(len(sources), -1, len(hours))
    # moon_position = source_elevation.reshape(len(sources), -1, len(hours))
    # source_elevation = source_elevation.reshape(len(sources), -1, len(hours))

    # print(source_elevation.reshape(len(sources), len(hours), -1))

    fig_1 = plt.figure()

    axis_1 = fig_1.add_subplot(111)

    for i, source in enumerate(sources):

        axis_1.imshow(source_elevation[i].T, aspect='auto', origin='lower',
                      label=source['name'])
        # plt.colorbar()

        axis_1.legend(loc='best')


    plt.show()


if __name__ == '__main__':
    main(location_filename='digicamscheduling/config/location_krakow.txt')
    # main()
