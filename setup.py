from setuptools import setup, find_packages


setup(
    name='digicamscheduling',
    version='0.1.0',
    packages=find_packages(),
    url='https://github.com/calispac/digicamscheduling',
    license='GNU GPL 3.0',
    author='Cyril Alispach',
    author_email='cyril.alispach@gmail.com',
    long_description=open('README.md').read(),
    description='A package for observation scheduling in gamma-ray astronomy',
    requires=['numpy', 'astropy', 'pyastronomy', 'matplotlib', 'scipy'],
    # bin=['bin/elevation.py'],
    entry_points={'console_scripts': [
        'digicamscheduling-elevation=digicamscheduling.bin.elevation:plot_elevation'],
    }
)