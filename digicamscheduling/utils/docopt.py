from schema import Schema, And, Use, Optional
from astropy.time import Time
import astropy.units as u

def convert_commandline_arguments(args):

    new_args = {}

    for key, val in args.items():
        new_args[key.replace('--', '')] = val


    schema = {'location_filename': Use(str),
              Optional('sources_filename'): Use(str),
              Optional('environment_filename'): Use(str),
              'start_date': Use(
                  lambda s: Time.now() if s is None else Time(s)),
              'end_date':
                  Use(lambda s: Time.now() + 1 * u.day if s is None else Time(
                      s)),
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

    print(args)


    return args
