from astroplan import Observer
from astroplan import FixedTarget

from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.time import Time

coordinates_krakow = {'lat': 50.090763 * u.deg, 'lon': 19.887956 * u.deg, 'height':  230 * u.m}

location = EarthLocation.from_geodetic(**coordinates_krakow)
telescope = Observer(location=location, name='SST-1M 1st prototype')

target_list = [FixedTarget.from_name('Crab')]

print(target_list[0])