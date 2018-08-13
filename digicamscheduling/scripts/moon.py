import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.core import moon
from digicamscheduling.io import reader
import digicamscheduling.display.plot as display
import matplotlib.pyplot as plt
from digicamscheduling.utils import time


if __name__ == '__main__':

    location_filename = 'digicamscheduling/config/' + 'location_krakow.txt'
    coordinates_krakow = reader.read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    start_date = Time('2018-07-27 12:00')
    end_date = Time('2018-07-28 12:00')
    time_steps = 1 * u.min
    date = time.compute_time(date_start=start_date, date_end=end_date,
                             time_steps=time_steps, location=location,
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


