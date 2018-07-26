import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.core import sun
from digicamscheduling.io import reader
import digicamscheduling.display.plot as display
import matplotlib.pyplot as plt
from digicamscheduling.utils import time


if __name__ == '__main__':

    location_filename = 'digicamscheduling/config/' + 'location_krakow.txt'
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    start_date = Time('2017-08-31 20:00')
    end_date = Time('2018-08-31 20:00')
    time_steps = 15 * u.min
    date = time.compute_time(date_start=start_date, date_end=end_date, time_steps=time_steps, location=location, only_night=False)
    sun_position = sun.compute_sun_position(date=date, location=location)

    fig_1 = plt.figure()
    fig_2 = plt.figure()
    fig_3 = plt.figure()
    axis_1 = fig_1.add_subplot(111)
    axis_2 = fig_2.add_subplot(111)
    axis_3 = fig_3.add_subplot(111, projection='polar')

    display.plot_azimuth(date, sun_position.az, axis=axis_1, label='Sun')
    display.plot_elevation(date, sun_position.alt, axis=axis_2, label='Sun')
    display.plot_trajectory(sun_position.az, sun_position.alt, axes=axis_3, label='Sun')

    plt.show()


