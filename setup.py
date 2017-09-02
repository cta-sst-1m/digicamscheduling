from distutils.core import setup


setup(
    name='digicamscheduling',
    version='0.1.0',
    packages=['digicamscheduling', 'digicamscheduling.core', 'digicamscheduling.display', 'digicamscheduling.io'],
    url='https://github.com/calispac/digicamscheduling',
    license='GNU GPL 3.0',
    author='Cyril Alispach',
    author_email='cyril.alispach@gmail.com',
    long_description=open('README.md').read(),
    #scripts=['bin/schedule.py', 'bin/visibility.py', 'bin/elevation.py', 'bin/crab_visibility.py', 'bin/source_position.py'],
    description='A package for observation scheduling in gamma-ray astronomy',
    requires=['numpy', 'astropy', 'pyastronomy', 'matplotlib', 'scipy'],
)