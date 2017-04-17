#!/usr/bin/env python

from setuptools import setup

with open('README.md') as file:
    long_description = file.read()

setup(name='pyquante2',
      version='0.1',
      description='Python Quantum Chemistry, version 2.0',
      long_description = long_description,
      author='Rick Muller',
      author_email='rpmuller@gmail.com',
      url='http://pyquante.sourceforge.net',
      install_requires=['numpy'],
      packages=['pyquante2',
                'pyquante2.basis',
                'pyquante2.dft',
                'pyquante2.geo',
                'pyquante2.graphics',
                'pyquante2.grid',
                'pyquante2.ints',
                'pyquante2.pt',
                'pyquante2.scf',
                'pyquante2.viewer',
                ],
      classifiers = ["Development Status :: 3 - Alpha",
                     "License :: OSI Approved :: BSD License",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering",
                     ],
      )
