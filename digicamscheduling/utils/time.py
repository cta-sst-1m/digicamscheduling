import numpy as np
import astropy.units as u
from digicamscheduling.core import sun


def compute_time(date_start, date_end, time_steps, location=None, only_night=True):

    duration = date_end.jd - date_start.jd
    time_bins = np.arange(0, duration, time_steps.to('day').value) * u.day
    time = date_start + time_bins

    if only_night:

        position_sun = sun.compute_sun_position(date=time, location=location)
        night = sun.compute_night(alt=position_sun.alt, type='astronomical')

        if (night == False).all():

            raise ValueError('The time specified does not'
                             ' contain any night period !')

        return time[night]

    else:

        return time


def compute_days(date_start, date_end):

    time_steps = 1 * u.day

    duration = date_end.jd - date_start.jd
    time_bins = np.arange(0, duration, time_steps.to('day').value) * u.day
    days = date_start + time_bins

    return days