#!/usr/bin/env python

import os
import glob

# We use numpy distutils to compile and wrap f90 code via f2py
import setuptools
    
with open('README.md', 'r') as fh:
    readme = fh.read()

with open('partycls/core/_version.py') as f:
    exec(f.read())

args = dict(name='partycls',
            version=__version__,
            description='Unsupervised learning of the structure of particulate systems',
            long_description=readme,
            long_description_content_type="text/markdown",
            author='Joris Paret',
            author_email='joris.paret@gmail.com',
            maintainer='Joris Paret',
            url='https://github.com/jorisparet/partycls',
            download_url='https://pypi.org/project/partycls',
            keywords=['physics', 'condensed matter', 'structure',
                      'machine learning', 'unsupervised', 'clustering',
                      'molecular dynamics', 'particle',
                      'descriptor', 'radial', 'angular', 'correlation'],
            packages=['partycls',
                      'partycls/core',
                      'partycls/descriptors'],
            install_requires=['numpy', 'sklearn'],
            license='GPLv3',
            setup_requires = ['numpy'],
            classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Development Status :: 5 - Production/Stable',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.6',                
                'Programming Language :: Python :: 3.7',
                'Programming Language :: Python :: 3.8',
                'Programming Language :: Python :: 3.9',
                'Topic :: Scientific/Engineering :: Physics'
                ]
)

try:
    from numpy.distutils.core import setup, Extension
    
    args["ext_modules"] = [Extension('partycls.descriptors.realspace_wrap',
                                     sources=['partycls/descriptors/realspace.f90'],
                                     extra_f90_compile_args=[]),
                           Extension('partycls.neighbors_wrap',
                                     sources=['partycls/neighbors.f90'],
                                     extra_f90_compile_args=[])]

except (ModuleNotFoundError, ImportError):
    from distutils.core import setup

setup(**args)
