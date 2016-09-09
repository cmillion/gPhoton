"""
.. module:: setup

   :synopsis: This script is used to setup the pip package.

.. moduleauthor:: Chase Million <chase.million@gmail.com>
"""

from __future__ import absolute_import
from os import path
# Always prefer setuptools over distutils
from setuptools import setup
# To use a consistent encoding
from gPhoton import __version__

HERE = path.abspath(path.dirname(__file__))

setup(
    name='gPhoton',
    version=__version__,
    description='The GALEX photon project.',
    url='https://github.com/cmillion/gPhoton',
    author='Chase Million, et al.',
    author_email='chase.million@gmail.com',
    zip_safe=False,
    license='AURA',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        ],
    keywords=['astronomy', 'galex', 'photometry', 'ultraviolet'],
    scripts=['bin/gPipeline', 'bin/gAperture', 'bin/gFind', 'bin/gMap'],
    packages=['gPhoton', 'gPhoton.cal'],
    install_requires=['numpy', 'scipy', 'requests>=2.4.0', 'pandas', 'astropy',],
    package_data={'gPhoton.cal' : ['cal/*.fits', 'cal/*.tbl']},
    include_package_data=True,
    download_url='https://archive.stsci.edu/prepds/gphoton/cal/gPhoton-{v}.tar'
    '.gz'.format(v=__version__),
)
