import time


def timeit(func):

    def timed(*args, **kw):

        ts = time.time()
        result = func(*args, **kw)
        te = time.time()
        print('[Time It] function: "{}({}, {})" \n run in {}s'.format(func.__name__, args, kw, te - ts))
        return result

    return timed
