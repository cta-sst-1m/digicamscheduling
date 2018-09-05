# digicam-extra
A repository containing shared config files for the various sst-1m projects

## catalog.json

JSON file containing the SST-1M source catalog. Coordinates are obtained from http://tevcat.uchicago.edu/ 

## catalog_stars.json

JSON file containing the SST-1M pointing stars. Coordinates are obtained from `astropy.coordinates.SkyCoord().from_name('NAME')`

## calibration__YYYYMMDD_*.yml

YAML file containing the calibration results
