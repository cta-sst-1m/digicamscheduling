import numpy as np
import astropy.units as u


def write_schedule(schedule, sources, dates, filename, units='deg'):

    observing = False
    previous_source_name = ''
    n_steps = len(dates)

    with open(filename, 'w') as file:

        for t, date in enumerate(dates):

            if t >= n_steps - 1:

                break

            date_string = '%s' % date.iso[0:-4]

            if np.any(schedule[..., t] == True):

                if not observing:

                    command_string = 'STARTUP'
                    file.write('%s' % (date - 90 * u.min).iso[0:-4] + '  ' +
                               command_string + '\n')
                    file.write('%s' % (date - 15 * u.min).iso[0:-4] + '  ' +
                               'READY' + '\n')

                    observing = True

                source_id = np.where(schedule[..., t] == True)[0][0]
                source = sources[source_id]
                if source['name'] != previous_source_name:

                    command_string = 'OBSERVING  ={"source": "%s", ' \
                                     '"dec": "%f", "ra": "%f"}'
                    command_string = command_string % (
                        source['name'], source['dec'].to(units).value,
                        source['ra'].to(units).value)
                    file.write(date_string + '  ' + command_string + '\n')

                    previous_source_name = source['name']

            else:

                if observing:
                    command_string = 'SHUTDOWN'
                    file.write(date_string + '  ' + command_string + '\n')
                    previous_source_name = ''

                observing = False

        # date += (time_bins[-1] - time_bins[-2])
        date_string = '%s' % date.iso[0:-4]
        command_string = 'SHUTDOWN'
        file.write(date_string + '  ' + command_string)
