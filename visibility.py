import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from digicamscheduling.core.scheduler import compute_source_visibility, find_priority_schedule, compute_schedule_efficiency, find_quality_schedule, find_dynamic_priority_quality_schedule
import digicamscheduling.core.sun as sun
import digicamscheduling.core.moon as moon
import digicamscheduling.core.gamma_source as gamma_source
import matplotlib.pyplot as plt
from digicamscheduling.io.reader import read_catalog, read_location


if __name__ == '__main__':

    sources_filename = 'digicamscheduling/config/' + 'catalog.txt'
    location_filename = 'digicamscheduling/config/' + 'location_krakow.txt'

    sources = read_catalog(sources_filename)
    coordinates_krakow = read_location(filename=location_filename)
    location = EarthLocation(**coordinates_krakow)

    time_bins = np.linspace(0, 1, num=24*6 + 1) * u.day
    time_interval = np.diff(time_bins)[0]
    start_date = Time('2017-10-30 12:00')

    sun_intensity = np.zeros(time_bins.shape)
    moon_intensity = np.zeros(time_bins.shape)
    source_intensity = np.zeros((len(sources), time_bins.shape[0]))

    for i, time in enumerate(time_bins):

        date = start_date + time
        sun_intensity[i] = sun.intensity(date=date, location=location)
        moon_intensity[i] = moon.intensity(date=date, location=location)

        for j, source in enumerate(sources):

            source_intensity[j, i] = gamma_source.intensity(date=date, location=location, ra=source['ra'], dec=source['dec'])

    sources_visibility = compute_source_visibility(source_intensity, sun_intensity, moon_intensity)
    sources_flux = np.array([source['flux'] for source in sources])
    t_eff = (np.sum(source_intensity, axis=1) * sources_flux * time_interval).to('hour')

    source_sorted = np.argsort(t_eff)
    fig, axis_1 = plt.subplots()

    for j in source_sorted:

        if t_eff[j] > 0:

            x = (start_date + time_bins).plot_date
            y = sources_visibility[j] * sources_flux[j]
            mask = (y >= 0)

            axis_1.plot_date(x[mask], y[mask], label=sources[j]['name'] + ' | $t_{eff}$' + ' = {0:0.02f}'.format(t_eff[j]), linestyle='-', marker='None')
            print(sources[j]['name'] + ' {0:0.02f}'.format(t_eff[j]))

    axis_1.set_xlabel('UTC time')
    axis_1.set_ylabel('effective flux [c.u.]')
    axis_1.set_ylim([0, 1.1])
    plt.legend(loc='best')
    plt.gcf().autofmt_xdate()
    plt.show()