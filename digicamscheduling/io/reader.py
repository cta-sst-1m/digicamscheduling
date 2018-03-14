import astropy.units as u
import pandas as pd


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


def read_ohp_weather_data(filename):
    df = pd.read_table(filename, header=1, skip_blank_lines=True, sep=',', converters={'Date': pd.to_datetime})
    names = df.columns.values.tolist()
    for i, name in enumerate(names):
        if name[0] == ' ':
            names[i] = name[1:]

    names = [name.replace(' ', '_') for name in names]
    names = [name.lower() for name in names]
    df.columns = names

    df = df.set_index('date')

    return df


if __name__ == '__main__':

    sources = read_catalog('../config/catalog.txt')

    # print(sources)

    read_ohp_weather_data('/data/datasets/CTA/weather_ohp/pluie_01012008_au_01012009.txt')