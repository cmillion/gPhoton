# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import glob

here = path.abspath(path.dirname(__file__))

(NAME, VERSION) = ('gPhoton', '1.23.4')
URL = 'https://github.com/cmillion/gPhoton'
AUTHOR = 'Chase Million, et al.'
EMAIL = 'chase.million@gmail.com'
LICENSE = 'AURA'

# Get the long description from the relevant file
#with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name = NAME,
    version = VERSION,
    description = 'The GALEX photon project.',
    #long_description = long_description,
    url = URL,
    author = AUTHOR,
    author_email = EMAIL,
    license = LICENSE,
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        ],
    keywords = 'astronomy galex photometry ultraviolet',
    packages = find_packages(exclude=['source.tests']),
    data_files=[('cal',[glob.glob('cal/*')]),
                ('e31000',[glob.glob('e31000/*')])],
    requires=['numpy','scipy','requests','pandas','astropy',],
)
