#!/usr/bin/env python

import os
import glob

# We use numpy distutils to compile and wrap f90 code via f2py
from numpy.distutils.core import setup, Extension

# Get the long description from README.md and try to convert it to
# reST. Adapted from https://bons.ai/blog/markdown-for-pypi
# This is necessary to have nicely formatted README on pypi until we
# understand how to make setuptools work nicely with f2py
try:
    from pypandoc import convert
    readme = convert('README.md', 'rst')
except (ImportError, OSError):
    try:
        readme = open('README.md', 'r').read()
    except:
        readme = ''

with open('partycls/core/_version.py') as f:
    exec(f.read())

args = dict(name='partycls',
            version=__version__,
            description='Structural clustering through correlation functions',
            long_description=readme,
            author='Joris Paret',
            author_email='joris.paret@umontpellier.fr',
            packages=['partycls',
                      'partycls/core',
                      'partycls/descriptor'],
            install_requires=['numpy', 'sklearn'],
            ext_modules=[Extension('partycls.descriptor.realspace_wrap',
                                   sources=['partycls/descriptor/realspace.f90'],
                                   extra_f90_compile_args=[])],
            license='GPLv3',
            setup_requires = ['numpy'],
            classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Development Status :: 4 - Beta',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.7',
                'Topic :: Scientific/Engineering :: Physics',
            ]
)

setup(**args)
