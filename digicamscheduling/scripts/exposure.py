from astropy.coordinates import EarthLocation, SkyCoord
from digicamscheduling.io import reader
from digicamscheduling.utils import time
from digicamscheduling.core import gamma_source
from digicamscheduling.display.plot import plot_sky_coord
from digicamscheduling.utils.exposure import compute_sidereal_exposure, compute_min_max_dec, normalize
from astropy.time import Time
import healpy as hp
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.linalg import norm
from astropy.table import Table
from matplotlib.cm import Reds, Reds_r

from astropy.utils.iers import conf
conf.auto_max_age = None


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
    x_ticks_label = to_ticks_label(np.flip(x_ticks))
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


def plot_labels_equatorial(axes=None, x_label='', y_label=''):
    axes.grid(True)
    x_ticks = np.linspace(-1, 1, num=13) * np.pi
    x_ticks_bis = np.copy(x_ticks)
    x_ticks_bis[0] = x_ticks_bis[0] - np.pi / 18
    x_ticks_bis[-1] = x_ticks_bis[-1] + np.pi / 12

    axes.set_xticks(x_ticks_bis)
    x_ticks = -(x_ticks - np.pi)
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


def matplotlib_to_healpy(theta, phi):

    theta = - theta + np.pi/2.
    phi = - phi + np.pi

    return theta, phi


def healpy_to_astropy(theta, phi):

    theta = - (theta - np.pi / 2.)

    return theta, phi

def healpy_to_matplotlib(theta, phi):

    theta = - (theta - np.pi / 2.)
    phi = - (phi - np.pi)

    return theta, phi


def astropy_to_matplotlib(theta, phi):

    x = - (phi - np.pi)

    return theta, x


def astropy_lat_lon_to_matplotlib(theta, phi):

    return theta, phi


def astropy_galactic_to_matplotlib(theta, phi):

    theta = np.atleast_1d(theta)
    phi = np.atleast_1d(phi)

    x_galactic = - phi + np.pi
    mask = x_galactic > 0
    x_galactic[mask] = x_galactic[mask] - np.pi
    x_galactic[~mask] = x_galactic[~mask] + np.pi
    phi = x_galactic

    return theta, phi

def astropy_galactic_to_matplotlib_bis(theta, phi):

    mask = phi > np.pi
    phi[mask] = -(phi[mask] - np.pi)
    return theta, -phi


if __name__ == '__main__':


    catalog = Table.read('gammacat.fits.gz')
    categories = np.array(np.unique(catalog['classes']))
    print(categories)
    count = np.arange(len(categories))

    for i, category in enumerate(categories):

        mask = catalog['classes'] == category
        count[i] = mask.sum()

    marker = ['+', '^', '*', 'x', 'o']
    colors = ['r', 'g', 'b', 'y', 'k']

    index = np.argsort(-count)
    n_classes = len(marker)

    sources_ra = np.array(catalog['ra'].to(u.rad).value)
    sources_dec = np.array(catalog['dec'].to(u.rad).value)
    sources_category = np.zeros(len(catalog['classes']))
    masked = np.zeros(len(catalog), dtype=bool)

    for i in range(n_classes):

        if i == n_classes - 1:
            sources_category[~masked] = i

            break

        mask = catalog['classes'] == categories[index[i]]
        masked += mask
        sources_category[mask] = i


    plt.figure()
    plt.bar(np.arange(len(count)), count)
    plt.xticks(np.arange(len(count)), categories)

    plt.figure()
    plt.hist(sources_category, bins=np.arange(8))

    plt.figure()
    plt.hist(sources_ra, bins='auto')
    plt.hist(sources_dec, bins='auto')

    ##  Plotting options
    projection = None
    projection = 'aitoff'
    # projection = 'mollweide'
    # projection = None
    alpha = 0.2
    inverse = 1

    fig_galactic = plt.figure()
    fig_equatorial = plt.figure()
    axes_galactic = fig_galactic.add_subplot(111, projection=projection)
    axes_equatorial = fig_equatorial.add_subplot(111, projection=projection)
    plot_labels(axes_galactic)
    plot_labels_equatorial(axes_equatorial)




    ##  Detector configuration
    lat_CTAS, lon_CTAS = -24.683428 * u.deg, -70.316344 * u.deg
    lat_CTAN, lon_CTAN = 28.756944 * u.deg, -17.886389 * u.deg
    lat_LHASSO, lon_LHASSO = 29.358611 * u.deg, 100.1375 * u.deg
    lat_HAWC, lon_HAWC = 18.995155 * u.deg, - 97.307153 * u.deg
    lat_north, lon_north = np.pi/2, 0
    lat_south, lon_south = -np.pi/2, 0

    fov_CTAS = 120 * u.deg
    fov_CTAN = 120 * u.deg
    fov_HAWC = 90 * u.deg
    fov_LHASSO = 90 * u.deg

    locations = []
    names = []
    fovs = []
    # location_CTAS = EarthLocation(lat=lat_CTAS, lon=lon_CTAS)
    # locations.append(location_CTAS)
    # names.append('CTA-South')
    # fovs.append(fov_CTAS)
    location_CTAN = EarthLocation(lat=lat_CTAN, lon=lon_CTAN)
    locations.append(location_CTAN)
    names.append('CTA-North')
    fovs.append(fov_CTAN)
    location_LHASSO = EarthLocation(lat=lat_LHASSO, lon=lon_LHASSO)
    locations.append(location_LHASSO)
    names.append('LHASSO')
    fovs.append(fov_LHASSO)
    # location_HAWC = EarthLocation(lat=lat_HAWC, lon=lon_HAWC)
    # locations.append(location_HAWC)
    # names.append('HAWC')
    # fovs.append(fov_HAWC)

    ## Healpix map and and coordinates
    n_pixel = 12 * 100 ** 2
    pixel = np.arange(n_pixel)
    n_side = hp.npix2nside(n_pixel)
    pixel_area = hp.nside2pixarea(n_side)

    dec, ra = hp.pix2ang(n_side, pixel)
    dir_sky = hp.ang2vec(dec, ra)
    y_equatorial, x_equatorial = healpy_to_matplotlib(dec, ra)
    dec_astropy, ra_astropy = healpy_to_astropy(dec, ra)
    sky = SkyCoord(ra=ra_astropy * u.rad, dec=dec_astropy * u.rad)
    coord = sky.transform_to('galactic')
    b = coord.b.radian
    l = coord.l.radian
    #x_equatorial = - (ra - np.pi)
    y_galactic, x_galactic = astropy_galactic_to_matplotlib(b, l)
    # x_galactic = np.unwrap(x_galactic)
    # _galactic[~mask] = np.pi - l[~mask]

    exposure = np.zeros(n_pixel)

    plt.figure()
    plt.hist(b, bins='auto', label='galactic latitude')
    plt.hist(dec, bins='auto', label='dec (healpy)')
    plt.hist(dec_astropy, bins='auto', label='dec (astropy)')
    plt.xlabel('[rad]')
    plt.legend(loc='best')

    plt.figure()
    plt.hist(l, bins='auto', label='galactic lattitude')
    plt.hist(ra, bins='auto', label='ra (healpy)')
    plt.hist(ra_astropy, bins='auto', label='ra (astropy)')
    plt.xlabel('[rad]')
    plt.legend(loc='best')

    label = ''

    for i in range(len(locations)):

        lat = locations[i].lat.to(u.rad).value
        max_zenith = (fovs[i]/2).to(u.rad).value
        exposure += compute_sidereal_exposure(dec_astropy, lat, max_zenith, normed=True)
        label += names[i] + ', '


    exposure /= len(locations) * np.degrees(pixel_area)
    mask = exposure > 0
    view = np.zeros(n_pixel)
    view[mask] = 1


    im = axes_equatorial.scatter(x_equatorial[mask], y_equatorial[mask],
                                 c=exposure[mask],
                                 label=label, cmap=Reds)


    # plt.colorbar(im, ax=axes_equatorial, label='Exposure [yr$^{-1}$]')
    # axes_equatorial.scatter(x_equatorial, dec, c=view, label=label, cmap=Reds)
    axes_galactic.scatter(x_galactic[mask], y_galactic[mask],
                          c=exposure[mask],
                          label=label, cmap=Reds)

    """
    plot_sky_coord(catalog['ra'].to(u.rad).value,
                   catalog['dec'].to(u.rad).value,
                   galactic=False,
                   color='k',
                   marker='.',
                   axes=axes_equatorial)
    plot_sky_coord(catalog['ra'].to(u.rad).value,
                   catalog['dec'].to(u.rad).value,
                   galactic=True,
                   color='k',
                   marker='.',
                   axes=axes_galactic)
                   
    """

    coord = SkyCoord(ra=lon_north * u.rad, dec=lat_north * u.rad)
    coord_galactic = coord.transform_to('galactic')
    x, y = astropy_galactic_to_matplotlib(coord_galactic.b.to(u.rad).value,
                                          coord_galactic.l.to(u.rad).value)
    axes_galactic.plot(y, x, marker='*', color='r',
                       linestyle='None', label='Earth North')


    coord = SkyCoord(ra=lon_south * u.rad, dec=lat_south * u.rad)
    coord_galactic = coord.transform_to('galactic')

    x, y = astropy_galactic_to_matplotlib(coord_galactic.b.to(u.rad).value,
                                          coord_galactic.l.to(u.rad).value)
    axes_galactic.plot(y, x, marker='*',
                       color='b',
                       linestyle='None', label='Earth South')

    axes_galactic.legend(loc='best')
    axes_equatorial.legend(loc='best')

    fig_galactic = plt.figure()
    fig_equatorial = plt.figure()
    axes_galactic = fig_galactic.add_subplot(111, projection=projection)
    axes_equatorial = fig_equatorial.add_subplot(111, projection=projection)
    plot_labels(axes_galactic)
    plot_labels_equatorial(axes_equatorial)

    """
    plot_sky_coord(catalog['ra'].to(u.rad).value,
                   catalog['dec'].to(u.rad).value,
                   galactic=False,
                   color='k',
                   marker='.',
                   axes=axes_equatorial)
    plot_sky_coord(catalog['ra'].to(u.rad).value,
                   catalog['dec'].to(u.rad).value,
                   galactic=True,
                   color='k',
                   marker='.',
                   axes=axes_galactic)
    """

    axes_equatorial.plot(lon_north, lat_north, marker='*', color='r',
                         linestyle='None', label='Earth North')
    axes_equatorial.plot(lon_south, lat_south, marker='*', color='b',
                         linestyle='None', label='Earth South')

    coord = SkyCoord(ra=lon_north * u.rad, dec=lat_north * u.rad)
    coord_galactic = coord.transform_to('galactic')
    x, y = astropy_galactic_to_matplotlib(coord_galactic.b.to(u.rad).value,
                                          coord_galactic.l.to(u.rad).value)
    axes_galactic.plot(y, x, marker='*', color='r',
                       linestyle='None', label='Earth North')
    coord = SkyCoord(ra=lon_south * u.rad, dec=lat_south * u.rad)
    coord_galactic = coord.transform_to('galactic')

    x, y = astropy_galactic_to_matplotlib(coord_galactic.b.to(u.rad).value,
                                          coord_galactic.l.to(u.rad).value)
    axes_galactic.plot(y, x, marker='*',
                       color='b',
                       linestyle='None', label='Earth South')

    exposure = np.zeros(n_pixel)
    view = np.zeros(n_pixel)

    for i in range(len(locations)):

        ra_obs = np.linspace(0,  2 * np.pi, 500)
        dec_obs = locations[i].lat.to(u.rad).value
        dec_obs = np.ones(ra_obs.shape) * dec_obs

        lon, lat = locations[i].lon.to(u.rad).value, locations[i].lat.to(u.rad).value
        max_zenith = (fovs[i]/2).to(u.rad).value

        coord_obs = SkyCoord(ra=lon * u.rad, dec=lat * u.rad)
        coord_obs_galactic = coord_obs.transform_to('galactic')
        lon_galactic = coord_obs_galactic.l.to(u.rad).value
        lat_galactic = coord_obs_galactic.b.to(u.rad).value

        coord = SkyCoord(ra=ra_obs * u.rad, dec=dec_obs * u.rad)
        coord_galactic = coord.transform_to('galactic')

        l_obs = coord_galactic.l.to(u.rad).value
        b_obs = coord_galactic.b.to(u.rad).value

        y, x = astropy_to_matplotlib(dec_obs, ra_obs)
        axes_equatorial.plot(x, y, label=names[i])
        y, x = astropy_lat_lon_to_matplotlib(lat, lon)
        axes_equatorial.plot(x, y, marker='o', color='k', linestyle='None')

        #x = b_obs
        #y = l_obs

        plt.figure()
        plt.hist(l_obs, bins='auto', label='l', alpha=0.3)
        plt.hist(b_obs, bins='auto', label='b', alpha=0.3)

        y, x = astropy_galactic_to_matplotlib(b_obs, l_obs)
        shift = len(x) - np.argmin(x) - 1
        x = np.roll(x, shift)
        y = np.roll(y, shift)

        plt.hist(x, bins='auto', label='l transformed', alpha=0.3)
        plt.hist(y, bins='auto', label='b transformed', alpha=0.3)
        plt.legend()

        axes_galactic.plot(x, y, label=names[i])
        y, x = astropy_lat_lon_to_matplotlib(lat_galactic, lon_galactic)
        axes_galactic.plot(-x, y, marker='o', color='k',linestyle='None')
        #axes_galactic.plot(np.pi/2, -np.pi/4, marker='o', color='r',linestyle='None')

        theta, phi = matplotlib_to_healpy(lat, lon)
        dir_detector = hp.ang2vec(theta, phi)
        angular_distance = (dir_sky * dir_detector).sum(axis=-1)

        angular_distance /= norm(dir_detector, axis=-1) * norm(dir_sky, axis=-1)
        angular_distance = np.arccos(angular_distance)
        # angular_distance = hp.rotator.angdist(dir_detector, dir_sky)
        mask = angular_distance <= max_zenith
        # temp = compute_sidereal_exposure(dec, lat, max_zenith)

        """

        for val in [np.pi / 3, np.pi/4, np.pi / 6]:

            mask = angular_distance <= val
            temp[mask] += 1
        """
        exposure = np.cos(angular_distance)
        exposure[~mask] = 0

        view += exposure


    mask = view > 0
    x = x_equatorial[mask]
    y = y_equatorial[mask]

    axes_equatorial.scatter(x, y, c=view[mask], cmap=Reds)
    x = x_galactic[mask]
    y = y_galactic[mask]
    axes_galactic.scatter(x, y, c=view[mask], cmap=Reds)

    axes_galactic.legend(loc='best')
    axes_equatorial.legend(loc='best')


    plt.show()