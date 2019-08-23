import numpy as np
import astropy.units as u


def compute_sidereal_exposure(dec, latitude, max_zenith, normed=False):

    lat = latitude
    max_z = max_zenith

    ks = np.cos(max_z) - np.sin(lat) * np.sin(dec)
    ks = ks / (np.cos(lat) * np.cos(dec))
    a = np.arccos(ks)
    a[ks > 1] = 0
    a[ks < -1] = np.pi

    exposure = np.cos(lat) * np.cos(dec) * np.sin(a)
    exposure = exposure + a * np.sin(lat) * np.sin(dec)

    if normed:
        exposure /= np.sum(exposure)

    return exposure


def compute_min_max_dec(latitude, max_zenith):

    dec_min = latitude - max_zenith
    dec_min = normalize(dec_min, lower=-np.pi/2., upper=np.pi/2.)

    dec_max = latitude + max_zenith
    dec_max = normalize(dec_max, lower=-np.pi/2., upper=np.pi/2.)

    return dec_min, dec_max


def normalize(num, lower=0.0, upper=360.0, b=False):
    """Normalize number to range [lower, upper) or [lower, upper].
    Parameters
    ----------
    num : float
        The number to be normalized.
    lower : float
        Lower limit of range. Default is 0.0.
    upper : float
        Upper limit of range. Default is 360.0.
    b : bool
        Type of normalization. See notes.
    Returns
    -------
    n : float
        A number in the range [lower, upper) or [lower, upper].
    Raises
    ------
    ValueError
      If lower >= upper.
    Notes
    -----
    If the keyword `b == False`, the default, then the normalization
    is done in the following way. Consider the numbers to be arranged
    in a circle, with the lower and upper marks sitting on top of each
    other. Moving past one limit, takes the number into the beginning
    of the other end. For example, if range is [0 - 360), then 361
    becomes 1. Negative numbers move from higher to lower
    numbers. So, -1 normalized to [0 - 360) becomes 359.
    If the keyword `b == True` then the given number is considered to
    "bounce" between the two limits. So, -91 normalized to [-90, 90],
    becomes -89, instead of 89. In this case the range is [lower,
    upper]. This code is based on the function `fmt_delta` of `TPM`.
    Range must be symmetric about 0 or lower == 0.
    Examples
    --------
    >>> normalize(-270,-180,180)
    90
    >>> import math
    >>> math.degrees(normalize(-2*math.pi,-math.pi,math.pi))
    0.0
    >>> normalize(181,-180,180)
    -179
    >>> normalize(-180,0,360)
    180
    >>> normalize(36,0,24)
    12
    >>> normalize(368.5,-180,180)
    8.5
    >>> normalize(-100, -90, 90, b=True)
    -80.0
    >>> normalize(100, -90, 90, b=True)
    80.0
    >>> normalize(181, -90, 90, b=True)
    -1.0
    >>> normalize(270, -90, 90, b=True)
    -90.0
    """
    from math import floor, ceil
    # abs(num + upper) and abs(num - lower) are needed, instead of
    # abs(num), since the lower and upper limits need not be 0. We need
    # to add half size of the range, so that the final result is lower +
    # <value> or upper - <value>, respectively.
    res = num
    if not b:
        if lower >= upper:
            raise ValueError("Invalid lower and upper limits: (%s, %s)" %
                             (lower, upper))

        res = num
        if num > upper or num == lower:
            num = lower + abs(num + upper) % (abs(lower) + abs(upper))
        if num < lower or num == upper:
            num = upper - abs(num - lower) % (abs(lower) + abs(upper))

        res = lower if res == upper else num
    else:
        total_length = abs(lower) + abs(upper)
        if num < -total_length:
            num += ceil(num / (-2 * total_length)) * 2 * total_length
        if num > total_length:
            num -= floor(num / (2 * total_length)) * 2 * total_length
        if num > upper:
            num = total_length - num
        if num < lower:
            num = -total_length - num

        res = num * 1.0  # Make all numbers float, to be consistent

    return res
