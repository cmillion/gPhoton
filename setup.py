"""
.. module:: setup

   :synopsis: This script is used to setup the pip package.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

# Always prefer setuptools over distutils
from setuptools import setup
# To use a consistent encoding
from os import path
from gPhoton import __version__

HERE = path.abspath(path.dirname(__file__))

setup(
    name='gPhoton',
    version=__version__,
    description='The GALEX photon project.',
    url='https://github.com/cmillion/gPhoton',
    author='Chase Million, et al.',
    author_email='chase.million@gmail.com',
    license='AURA',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        ],
    keywords=['astronomy', 'galex', 'photometry', 'ultraviolet'],
    scripts=['bin/gPipeline', 'bin/gAperture', 'bin/gFind', 'bin/gMap'],
    packages=['gPhoton', 'gPhoton.cal'],
    install_requires=['numpy>=2.4.0', 'scipy', 'requests', 'pandas', 'astropy',],
    package_data={'gPhoton.cal' : ['cal/*.fits', 'cal/*.tbl']},
    include_package_data=True,
    download_url='https://archive.stsci.edu/prepds/gphoton/cal/gPhoton-{v}.tar'
    '.gz'.format(v=__version__),
)
