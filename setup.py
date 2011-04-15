#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='pyquante2',
      version='0.1',
      description='Python Quantum Chemistry, version 2.0',
      author='Rick Muller',
      author_email='rpmuller@gmail.com',
      url='http://pyquante.sourceforge.net',
      install_requires=['numpy>=1.3'],
      packages=find_packages(),
      test_suite='pyquante2.tests.test_all',
      )
