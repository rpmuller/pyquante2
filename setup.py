#!/usr/bin/env python

from distutils.core import setup

setup(name='pyquante2',
      version='0.1',
      description='Python Quantum Chemistry, version 2.0',
      author='Rick Muller',
      author_email='rpmuller@gmail.com',
      url='http://pyquante.sourceforge.net',
      #install_requires=["numpy>=1.1"], # only works with setuptools, I guess.
      packages=['pyquante2',
                'pyquante2.basis',
                'pyquante2.ints',
                ],
      )
