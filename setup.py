from setuptools import setup, find_packages

command_line_applications = {'console_scripts': [
    'digicamscheduling-catalog=digicamscheduling.scripts.catalog:entry',
    'digicamscheduling-elevation=digicamscheduling.scripts.elevation:entry',
    'digicamscheduling-moon=digicamscheduling.scripts.moon:entry',
    'digicamscheduling-observability=digicamscheduling.scripts.observability:entry',
    'digicamscheduling-sun=digicamscheduling.scripts.sun:entry',
    'digicamscheduling-schedule=digicamscheduling.scripts.schedule:entry',

],
       }

setup(
    name='digicamscheduling',
    version='0.2.1',
    packages=find_packages(),
    url='https://github.com/calispac/digicamscheduling',
    license='GNU GPL 3.0',
    author='Cyril Alispach',
    author_email='cyril.alispach@gmail.com',
    long_description=open('README.md').read(),
    package_data={'': ['config/*']},
    description='A package for observation scheduling in gamma-ray astronomy',
    requires=['numpy', 'astropy', 'pyastronomy', 'matplotlib', 'scipy',
              'docopt'],
    entry_points=command_line_applications,
)
