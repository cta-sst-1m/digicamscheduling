from digicamscheduling.io import reader
from digicamscheduling.core import environement

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


def plot_trees(azimuth, elevation, axis=None, **kwargs):

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111, projection='polar')

    r = azimuth.to('radian').value
    theta = (90 * u.deg - elevation).to('deg').value

    axis.fill_between(r, 90, theta, **kwargs)
    axis.set_rmax(90)
    axis.set_theta_zero_location("N")
    axis.set_theta_direction(-1)
    axis.set_yticks(np.arange(0, 90 + 10, 10))
    axis.set_yticklabels(axis.get_yticks()[::-1])
    axis.legend()

    return axis


if __name__ == '__main__':

    filename = 'digicamscheduling/config/environmental_limitation.txt'
    alt, az = reader.read_environmental_limits(filename) * u.deg
    env_limits_function = environement.interpolate_environmental_limits(alt,
                                                                        az)
    az = np.linspace(0, 360, num=1000) * u.deg
    alt = env_limits_function(az) * u.deg

    plot_trees(az, alt, facecolor='black', alpha=0.8)
    plt.show()
