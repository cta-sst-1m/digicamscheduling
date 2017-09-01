import astropy.units as u


def read_catalog(filename):
    sources = []

    with open(filename) as file:
        for line in file:
            line = line.split('  ')
            sources.append({'name': line[0], 'ra': float(line[1]) * u.deg, 'dec': float(line[2]) * u.deg, 'flux': float(line[3].rstrip())})

    return sources


def read_location(filename):

    location = {}

    with open(filename) as file:
        for line in file:
            line = line.split('  ')
            location['lat'] = float(line[0]) * u.deg
            location['lon'] = float(line[1]) * u.deg
            location['height'] = float(line[2].rstrip()) * u.m

    return location


if __name__ == '__main__':


    sources = read_catalog('../config/fact_catalog.txt')

    print(sources)