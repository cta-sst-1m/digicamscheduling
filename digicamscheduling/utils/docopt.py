from schema import Schema, And, Use, Optional
from astropy.time import Time
import astropy.units as u
import os
from pkg_resources import resource_filename

PACKAGE_NAME = 'digicamscheduling'
CONFIG_FOLDER = 'config/digicamextra/'

DEFAULT_LOCATION_FILENAME = os.path.join(CONFIG_FOLDER, 'location_krakow.txt')
DEFAULT_LOCATION_FILENAME = resource_filename(PACKAGE_NAME,
                                              DEFAULT_LOCATION_FILENAME)

DEFAULT_SOURCE_FILENAME = os.path.join(CONFIG_FOLDER, 'catalog.json')
DEFAULT_SOURCE_FILENAME = resource_filename(PACKAGE_NAME,
                                            DEFAULT_SOURCE_FILENAME)

DEFAULT_ENVIRONMENT_FILENAME = os.path.join(CONFIG_FOLDER,
                                            'environmental_limitation.txt')
DEFAULT_ENVIRONMENT_FILENAME = resource_filename(PACKAGE_NAME,
                                                 DEFAULT_ENVIRONMENT_FILENAME)


def convert_commandline_arguments(args):

    new_args = {}

    for key, val in args.items():
        new_args[key.replace('--', '')] = val

    schema = {'location_filename': Use(get_location_filename),
              Optional('sources_filename'): Use(get_sources_filename),
              Optional('environment_filename'): Use(get_environment_filename),
              'start_date': Use(compute_start_date),
              'end_date': Use(compute_end_date),
              'time_step': And(Use(float), Use(lambda t: t * u.minute)),
              'output_path': Use(lambda s: s if not s else str(s)),
              'help': Use(bool),
              Optional('hide'): Use(bool),
              Optional('threshold'): Use(float),
              Optional('use_moon'): Use(bool),
              }

    schema = Schema(schema)
    args = schema.validate(new_args)

    del args['help']

    return args


def compute_start_date(text):

    if text is None:

        return Time.now()

    elif isinstance(text, str):

        return Time(text)

    else:

        raise TypeError


def compute_end_date(text):
    if text is None:

        return Time.now() + 24 * u.hour

    elif isinstance(text, str):

        return Time(text)

    else:

        raise TypeError


def get_sources_filename(text):

    if text is None:

        return DEFAULT_SOURCE_FILENAME

    else:

        return str(text)


def get_environment_filename(text):

    if text is None:

        return DEFAULT_ENVIRONMENT_FILENAME

    else:

        return str(text)


def get_location_filename(text):

    if text is None:

        return DEFAULT_LOCATION_FILENAME

    else:

        return str(text)
