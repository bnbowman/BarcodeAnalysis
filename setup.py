import os
import sys

from setuptools import setup, find_packages

__author__ = 'Brett Bowma'

# Check that a valid, modern version of Python is being used
if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    print "pbhml requires Python 2.7"
    sys.exit(-1)

# Read the version information
vFile = 'pbbarcode2/_version.py'
if os.path.exists(vFile):
    lines = open(vFile, 'r').read().splitlines()
    for line in lines:
        elts = line.split('=')
        elts = [e.strip() for e in elts]
        if len(elts) == 2 and elts[0] == '__version__':
            _ReadVersion = elts[1].replace('\'', '').replace('\"', '')
            break
else:
    _ReadVersion = '0.0.0'

# Run Setup
setup(
    name = 'pbbarcode2',
    version=_ReadVersion,
    author='Brett Bowman',
    author_email='bbowman@pacificbiosciences.com',
    license='LICENSE.txt',
    packages = find_packages('pbhml'),
    package_dir = {'':'pbhml'},
    zip_safe = False,
    install_requires=[
        'numpy >= 1.6.0',
        'h5py >= 1.3.0',
        'ConsensusCore >= 0.8.8'
        ]
    )