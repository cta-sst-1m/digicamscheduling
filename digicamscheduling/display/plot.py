import matplotlib.pyplot as plt


def plot_azimuth(date, azimut, **kwargs):

    fig = plt.figure()
    axis = fig.add_subplot(111)

    axis.plot_date(date.plot_date, azimut, **kwargs)
    axis.set_xlabel('UTC time')
    axis.set_ylabel('azimuth [deg]')
    axis.legend(loc='upper left')
    plt.gcf().autofmt_xdate()

    return axis


def plot_elevation(date, elevation, **kwargs):

    fig = plt.figure()
    axis = fig.add_subplot(111)

    axis.plot_date(date.plot_date, elevation, **kwargs)
    axis.set_xlabel('UTC time')
    axis.set_ylabel('elevation [deg]')
    axis.legend(loc='upper left')
    plt.gcf().autofmt_xdate()

    return axis
