import numpy as np
import astropy.units as u
from digicamscheduling.core import sun


def compute_time(date_start, date_end, time_steps, location, only_night=True):

    duration = date_end.jd - date_start.jd
    time_bins = np.arange(0, duration, time_steps.to('day').value) * u.day
    time = date_start + time_bins

    if only_night:

        position_sun = sun.compute_sun_position(date=time, location=location)
        night = sun.compute_night(alt=position_sun.alt, type='astronomical')

        night_bins = np.arange(night.shape[0])[night]

        for index, night_bin in enumerate(night_bins):

            if night_bins[min(index+1, len(night_bins) -1)] - night_bin > 1:

                night[night_bin + 1] = True

        return time[night]

    else:

        return time
