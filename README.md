# digicamscheduling

Scheduling package for SST-1M observations

## Getting Started

### Prerequisites

```
Numpy, Scipy, Astropy, PyAstronomy, Matplolib, Pandas, Tqdm, Cython, Docopt
```

### Installing (with Anaconda)

```
git clone https://github.com/cta-sst-1m/digicamscheduling
cd digicamscheduling
conda env create -f environment.yml
source activate digicamscheduling
python setup.py install
```
Try one of command line scripts

```
digicamscheduling-catalog
digicamscheduling-observability
digicamscheduling-elevation
```

Use the option `--help` to see how to run the scripts

## Considerations

Please consider the time and time steps you use ! Memory is limited !!!
so you cannot compute for the entire catalog for the entire year
with 1 seconde time step.

## Know issues

Most of the issues appear when:

1. Date does not exists : e.g. 2018-02-31 
2. Time step is too low --> Memory error
3. Period is shorter than a day
4. Period does not include night time (Sun elevation < -12 deg)
