import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np


def plot_azimuth(date, azimuth, axis=None, **kwargs):

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)

    axis.plot_date(date.plot_date, azimuth, linestyle='None', marker='.', **kwargs)
    axis.set_xlabel('UTC time')
    axis.set_ylabel('azimuth [deg]')
    axis.legend(loc='best')
    plt.gcf().autofmt_xdate()
    for tick in axis.get_xticklabels():
        tick.set_rotation(45)
    return axis


def plot_elevation(date, elevation, axis=None, linestyle='None', marker='.', **kwargs):

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)

    axis.plot_date(date.plot_date, elevation, linestyle=linestyle, marker=marker, **kwargs)
    axis.set_xlabel('UTC time')
    axis.set_ylabel('elevation [deg]')
    axis.set_ylim([0, 90])
    axis.legend(loc='best')
    plt.gcf().autofmt_xdate()
    for tick in axis.get_xticklabels():
        tick.set_rotation(45)

    return axis


def plot_trajectory(azimuth, elevation, axis=None, **kwargs):

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111, projection='polar')

    axis.plot(azimuth.to('radian'), (90 * u.deg - elevation).to('deg'), linestyle='None', marker='.', **kwargs)
    axis.set_rmax(90)
    axis.set_theta_zero_location("N")
    axis.set_theta_direction(-1)
    axis.set_yticks(np.arange(0, 90 + 10, 10))
    axis.set_yticklabels(axis.get_yticks()[::-1])
    # axis.set_xticks([0, 45, 90, 135, 180, 225, 270, 315])
    # axis.set_xticklabels([0, 45, 90, 135, 180, -135, -90, -45])
    axis.legend()

    return axis


def plot_moon_phase(date, phase, axis=None, **kwargs):

    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)

    axis.plot_date(date.plot_date, phase, linestyle='None', marker='.', **kwargs)
    axis.set_xlabel('UTC time')
    axis.set_ylabel('phase [a.u.]')
    axis.set_ylim([0, 1.1])
    axis.legend(loc='best')
    plt.gcf().autofmt_xdate()
    for tick in axis.get_xticklabels():
        tick.set_rotation(45)

    return axis
