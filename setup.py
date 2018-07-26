from setuptools import setup, find_packages

command_line_applications = {'console_scripts':
                                 ['digicamscheduling-test = digicamscheduling.bin.test:main',
                                  ],
                             'gui_scripts':
                                 ['digicamscheduling-elevation = digicamscheduling.bin.elevation:main']
                             }
setup(
    name='digicamscheduling',
    # version='0.2.1',
    packages=find_packages(),
    url='https://github.com/calispac/digicamscheduling',
    license='GNU GPL 3.0',
    author='Cyril Alispach',
    author_email='cyril.alispach@gmail.com',
    long_description=open('README.md').read(),
    description='A package for observation scheduling in gamma-ray astronomy',
    requires=['numpy', 'astropy', 'pyastronomy', 'matplotlib', 'scipy'],
    entry_points=command_line_applications,
)
