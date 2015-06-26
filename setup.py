# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import glob
from gPhoton import __version__

here = path.abspath(path.dirname(__file__))

setup(
    name = 'gPhoton',
    version = __version__,
    description = 'The GALEX photon project.',
    url = 'https://github.com/cmillion/gPhoton',
    author = 'Chase Million, et al.',
    author_email = 'chase.million@gmail.com',
    license = 'AURA',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        ],
    keywords = ['astronomy', 'galex', 'photometry', 'ultraviolet'],
    scripts=['bin/gPipeline','bin/gAperture','bin/gFind','bin/gMap'],
    packages=['gPhoton','gPhoton.cal'],
    install_requires=['numpy','scipy','requests','pandas','astropy',],
    zip_safe = True,
    include_package_data = True,
    download_url = 'http://github.com/cmillion/gPhoton/tarball/#,gPhoton-{v}'.format(v=__version__),
)
