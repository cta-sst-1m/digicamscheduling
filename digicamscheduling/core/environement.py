import astropy.units as u


def compute_forbiden_zone(alt, az):

    return alt > 45 * u.deg