import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
import healpy as hp
from astropy.coordinates import SkyCoord


def plot_source(date, position, y_label='', axes=None, ylim=None, **kwargs):
    """
    fig = axes.get_figure()
    temp = date.mjd - date.mjd.min()
    n_days = int(temp.max())

    fig, axes = plt.subplots(1, n_days, sharey=True)
    axes[0].set_xlabel('UTC time')
    axes[0].set_ylim(ylim)

    for i, day in enumerate(range(1, n_days+1, 1)):

        t = date[temp <= day].plot_date
        y = position[temp <= day]

        axes[i].plot_date(t, y, **kwargs)
        axes[i].set_xlim(t.min(), t.max())
        axes[i].set_ylabel(y_label)
        if ylim is not None:
        for tick in axes[i].get_xticklabels():
            tick.set_rotation(45)
    axes[-1].legend(loc='best')
    fig.autofmt_xdate()

    """
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    try:

        y = position.value

    except Exception:

        y = position

        pass

    if ylim is not None:
        axes.set_ylim(ylim)
        mask = (y > ylim[0]) * (y < ylim[1])
        y = np.ma.masked_array(y, ~mask)

    axes.plot_date(date.plot_date, y, **kwargs)
    axes.set_xlabel('UTC time')
    axes.set_ylabel(y_label)

    plt.gcf().autofmt_xdate()
    for tick in axes.get_xticklabels():
        tick.set_rotation(45)

    return axes

    """
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(111)

    title = '\nTime start : {} Time end : {}'.format(date.min(), date.max())

    axes.set_title(title, fontsize=20)

    t = np.linspace(0, 1, num=len(date))
    axes.plot(t, position, **kwargs)
    axes.set_xlabel('time []')
    axes.set_ylabel(y_label)
    if ylim is not None:
        axes.set_ylim(ylim)
    axes.legend(loc='best')
    # plt.gcf().autofmt_xdate()
    # for tick in axes.get_xticklabels():
    #    tick.set_rotation(45)

    return axes
    




    if axes is None:

        fig = plt.figure()
        axes = fig.add_subplot(111)

    else:

        fig = axes.get_figure()

    temp = date.mjd - date.mjd.min()
    n_days = int(temp.max())

    days = np.arange(0, n_days, 1)

    print(temp)
    mask = np.searchsorted(temp, days)
    mask = date[mask].plot_date
    print(mask, temp)
    print(mask.shape, temp)

    xlims = [(mask[i], mask[i+1]) for i in range(n_days-1)]
    xlims = tuple(xlims)
    print(xlims)

    baxes = brokenaxes(xlims=xlims, ylims=None, d=.015, tilt=45,
                       subplot_spec=None, fig=fig, despine=True,)

    baxes.plot_date(date.plot_date, position, **kwargs)
    baxes.set_xlabel('UTC time')
    baxes.set_ylabel(y_label)
    if ylim is not None:
        axes.set_ylim(ylim)
    baxes.legend(loc='best')
    fig.autofmt_xdate()
    for tick in axes.get_xticklabels():
        tick.set_rotation(45)

    return baxes
    
    """


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


def to_ticks_label(x):
    y = []
    x = np.degrees(x)

    for i, val in enumerate(x):

        if (i == 0 or i == len(x) - 1):

            y.append(str(int(val)) + '°')
        else:
            y.append('')

    return y


def plot_labels(axes=None, x_label='', y_label=''):
    axes.grid(True)
    x_ticks = np.linspace(1, -1, num=13) * np.pi
    x_ticks_bis = np.copy(x_ticks)
    x_ticks_bis[0] = x_ticks_bis[0] + np.pi / 12
    x_ticks_bis[-1] = x_ticks_bis[-1] - np.pi / 18

    axes.set_xticks(x_ticks_bis)
    # x_ticks = -(x_ticks - np.pi)
    x_ticks_label = to_ticks_label(x_ticks)
    axes.set_xticklabels(x_ticks_label, fontdict={'verticalalignment': 'center',
                                                  'horizontalalignment': 'center'}
                         )

    y_ticks = np.linspace(-1 / 2, 1 / 2, num=7) * np.pi
    y_ticks_label = to_ticks_label(y_ticks)
    axes.set_yticks(y_ticks)
    axes.set_yticklabels(y_ticks_label, fontdict={'verticalalignment': 'center',
                                                  'horizontalalignment': 'center'})

    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)

    return axes


def astropy_galactic_to_matplotlib(theta, phi):

    theta = np.atleast_1d(theta)
    phi = np.atleast_1d(phi)

    x_galactic = - phi + np.pi
    mask = x_galactic > 0
    x_galactic[mask] = x_galactic[mask] - np.pi
    x_galactic[~mask] = x_galactic[~mask] + np.pi
    phi = x_galactic

    return theta, phi


def astropy_to_matplotlib(theta, phi):

    y = np.atleast_1d(theta)
    x = np.atleast_1d(phi)
    x = - x + np.pi

    return y, x


def plot_sky_coord(ra, dec, galactic=False, axes=None, **kwargs):

    if axes is None:

        fig = plt.figure()
        axes = fig.add_subplot(111)

    y, x = astropy_to_matplotlib(dec, ra)

    if galactic:

        sky = SkyCoord(ra=ra * u.rad, dec=dec * u.rad, frame='icrs')
        coord = sky.transform_to('galactic')
        b = coord.b.radian
        l = coord.l.radian
        y, x = astropy_galactic_to_matplotlib(b, l)

    axes.scatter(x, y, **kwargs)

    return axes