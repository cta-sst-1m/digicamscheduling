from schema import Schema, And, Use
from astropy.time import Time
import astropy.units as u

def convert_commandline_arguments(args):

    schema = {'--location_filename': Use(str),
              '--sources_filename': Use(str),
              '--environment_filename': Use(str),
              '--start_date': Use(
                  lambda s: Time.now() if s is None else Time(s)),
              '--end_date':
                  Use(lambda s: Time.now() + 1 * u.day if s is None else Time(
                      s)),
              '--time_step': And(Use(float), Use(lambda t: t * u.minute)),
              '--output_path': Use(lambda s: s if not s else str(s)),
              '--help': Use(bool),
              '--show': Use(bool),
              '--threshold': Use(float),
              '--use_moon': Use(bool),
              }

    schema = Schema(schema)
    args = schema.validate(args)

    for key, val in args.items():

        args[key.strip('--')] = args.pop(key)

    del args['help']

    return args
