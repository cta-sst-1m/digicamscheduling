from distutils.core import setup

setup(
    name='digicam_scheduling',
    version='0.1.0',
    packages=['digicam_scheduling',],
    url='https://github.com/calispac/digicam_scheduling',
    license='GNU GPL 3.0',
    author='Cyril Alispach',
    author_email='cyril.alispach@gmail.com',
    scripts=['digicam_scheduling/schedule.py', 'digicam_scheduling/visibility.py', 'digicam_scheduling/elevation.py', 'digicam_scheduling/crab_visibility.py'],
    description='A package for observation scheduling in gamma-ray astronomy',
    requires=['numpy', 'astropy', 'PyAstronomy', 'matplotlib', 'scipy'],
)