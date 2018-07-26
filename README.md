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

Please consider the time and time steps you use ! Memory is limited !!!
so you cannot compute for the entire catalog for the entire year
with 1 seconde time step.