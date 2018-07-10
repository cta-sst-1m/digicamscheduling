from digicamscheduling.io import reader
from digicamscheduling.core import environement
from digicamscheduling.display import plot

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


if __name__ == '__main__':

    filename = 'digicamscheduling/config/environmental_limitation.txt'
    alt, az = reader.read_environmental_limits(filename)
    env_limits_function = environement.interpolate_environmental_limits(alt,
                                                                        az)

    x = np.linspace(0, 360, num=1000) * u.deg

    print(env_limits_function(x))

    plot.plot_trajectory(x, env_limits_function(x) * u.deg)
    plot.plot_trajectory(az * u.deg, alt * u.deg)
    plt.show()
