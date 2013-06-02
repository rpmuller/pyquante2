#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = {}
ext_modules = []

if use_cython:
    ext_modules += [
        Extension("pyquante2.cutils",["cython/cutils.pyx"])
        ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("pyquante2.cutils",["cython/cutils.c"])
        ]

setup(name='pyquante2',
      version='0.1',
      description='Python Quantum Chemistry, version 2.0',
      author='Rick Muller',
      author_email='rpmuller@gmail.com',
      url='http://pyquante.sourceforge.net',
      packages=['pyquante2',
                'pyquante2.basis',
                'pyquante2.geo',
                'pyquante2.ints',
                ],
      cmdclass = cmdclass,
      ext_modules = ext_modules,
      )
