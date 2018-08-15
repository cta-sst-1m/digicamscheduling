import astropy.units as u
import pandas as pd
import numpy as np
import json


def read_catalog(filename):

    with open(filename, 'r') as file:

        sources = json.load(file)['sources']

    for source in sources:

        source['ra'] = source['ra'] * u.deg
        source['dec'] = source['dec'] * u.deg
        source['weight'] = float(source['weight'])

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
    df = pd.read_table(filename, header=1, skip_blank_lines=True, sep=',',
                       converters={'Date': pd.to_datetime})
    names = df.columns.values.tolist()
    for i, name in enumerate(names):
        if name[0] == ' ':
            names[i] = name[1:]

    names = [name.replace(' ', '_') for name in names]
    names = [name.lower() for name in names]
    df.columns = names

    df = df.set_index('date')

    return df


def read_environmental_limits(filename):

    data = np.loadtxt(filename)
    az = data[:, 0]
    alt = data[:, 1]

    return alt, az


if __name__ == '__main__':

    sources_2 = read_catalog('../config/catalog.json')

    read_ohp_weather_data('/data/datasets/CTA/weather_ohp/'
                          'pluie_01012008_au_01012009.txt')
