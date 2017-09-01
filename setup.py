from setuptools import setup, find_packages

setup(
    name='digicam_scheduling',
    version='0.1.0',
    packages=['inout', 'core', 'display', ''],
    url='https://github.com/calispac/digicam_scheduling',
    license='GNU GPL 3.0',
    author='Cyril Alispach',
    author_email='cyril.alispach@gmail.com',
    scripts=['schedule.py', 'visibility.py', 'elevation.py', 'crab_visibility.py'],
    description='A package for observation scheduling in gamma-ray astronomy',
    requires=['numpy', 'astropy', 'PyAstronomy', 'matplotlib', 'scipy']
)