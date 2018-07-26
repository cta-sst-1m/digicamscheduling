import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from matplotlib.dates import DateFormatter


def plot_source(date, position, y_label='', axes=None, ylim=None, **kwargs):

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    axes.plot_date(date.plot_date, position, **kwargs)
    axes.set_xlabel('UTC time')
    axes.set_ylabel(y_label)
    if ylim is not None:
        axes.set_ylim(ylim)
    axes.legend(loc='best')
    plt.gcf().autofmt_xdate()
    for tick in axes.get_xticklabels():
        tick.set_rotation(45)

    return axes


def plot_azimuth(date, azimuth, axes=None, **kwargs):

    return plot_source(date, azimuth, y_label='azimuth [deg]', axes=axes,
                       **kwargs)


def plot_elevation(date, elevation, axes=None, **kwargs):

    return plot_source(date, elevation, y_label='elevation [deg]', axes=axes,
                       ylim=[0, 90], linestyle='None', marker='.', **kwargs)


def plot_trajectory(azimuth, elevation, axes=None, **kwargs):

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111, projection='polar')

    axes.plot(azimuth.to('radian'), (90 * u.deg - elevation).to('deg'),
              linestyle='None', marker='.', **kwargs)
    axes.set_rmax(90)
    axes.set_theta_zero_location("N")
    axes.set_theta_direction(-1)
    axes.set_yticks(np.arange(0, 90 + 10, 10))
    axes.set_yticklabels(axes.get_yticks()[::-1])
    # axes.set_xticks([0, 45, 90, 135, 180, 225, 270, 315])
    # axes.set_xticklabels([0, 45, 90, 135, 180, -135, -90, -45])
    axes.legend()

    return axes


def plot_moon_phase(date, phase, axes=None, **kwargs):

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    axes.plot_date(date.plot_date, phase, linestyle='None', marker='.',
                   **kwargs)
    axes.set_xlabel('UTC time')
    axes.set_ylabel('phase [a.u.]')
    axes.set_ylim([0, 1.1])
    axes.legend(loc='best')
    plt.gcf().autofmt_xdate()
    for tick in axes.get_xticklabels():
        tick.set_rotation(45)

    return axes


def plot_source_2d(source_elevation, coordinates, source=None,
                   c_label='elevation [deg]', extent=None,
                   cmap=plt.get_cmap('RdYlGn'), axes=None, **kwargs):

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    cmap = cmap
    cmap.set_bad(color='k', alpha=1.)

    title = 'Site : ({}, {}, {}) [lat, lon, alt] deg'.format(
        coordinates['lat'], coordinates['lon'], coordinates['height']
    )

    if source is not None:

        title += '\nSource : {} ({}, {}) [ra, dec] deg'.format(source['name'],
                                                               source['ra'],
                                                                source['dec'])

    axes.set_title(title, fontsize=20)

    ax = axes.imshow(source_elevation.T, aspect='auto', origin='lower',
                     extent=extent, cmap=cmap, **kwargs)

    axes.xaxis_date()
    axes.xaxis.set_minor_formatter(DateFormatter('%Y-%m-%d'))
    axes.set_xlabel('date')
    axes.set_ylabel('hour [UTC]')
    axes.set_yticks(np.arange(0, 24, 1))

    fig = axes.get_figure()
    fig.autofmt_xdate()
    fig.colorbar(ax, label=c_label, extend='both')

    return axes


def plot_sun_2d(sun_elevation, coordinates, extent, axes=None, **kwargs):

    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    axes = plot_source_2d(sun_elevation, coordinates, extent=extent,
                               vmin=-90, vmax=90, axes=axes,
                               cmap=plt.get_cmap('RdYlGn_r'),
                               c_label='Sun elevation [deg]', **kwargs)

    cs = axes.contour(sun_elevation.T, levels=[-18., -12., -6., -0.],
                      extent=extent, cmap='binary_r')

    contour_labels = ['Astronomical', 'Nautical', 'Civil',
                      'Horizon']
    fmt = dict(zip(cs.levels, contour_labels))
    axes.clabel(cs, fmt=fmt)

    return axes
