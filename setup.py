#!/usr/bin/env python

from setuptools import setup
from setuptools.extension import Extension
import numpy as np

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
        Extension("pyquante2.cints.one", ["cython/cone_wrap.pyx", "cython/cints.c"]),
        Extension("pyquante2.cints.two", ["cython/ctwo_wrap.pyx", "cython/cints.c"]),
        Extension(
            "pyquante2.cints.hgp",
            ["cython/chgp_wrap.pyx", "cython/chgp.c", "cython/cints.c"],
        ),
        Extension("pyquante2.cints.rys", ["cython/crys_wrap.pyx", "cython/crys.c"]),
        Extension(
            "pyquante2.cbecke", ["cython/cbecke.pyx"], include_dirs=[np.get_include()]
        ),
    ]
    cmdclass.update({"build_ext": build_ext})
else:
    ext_modules += [
        Extension("pyquante2.cints.one", ["cython/cone_wrap.c", "cython/cints.c"]),
        Extension("pyquante2.cints.two", ["cython/ctwo_wrap.c", "cython/cints.c"]),
        Extension(
            "pyquante2.cints.hgp",
            ["cython/chgp_wrap.c", "cython/chgp.c", "cython/cints.c"],
        ),
        Extension("pyquante2.cints.rys", ["cython/crys_wrap.c", "cython/crys.c"]),
    ]

with open("README.md") as file:
    long_description = file.read()

setup(
    name="pyquante2",
    install_requires=["numpy"],
    packages=[
        "pyquante2",
        "pyquante2.basis",
        "pyquante2.dft",
        "pyquante2.geo",
        "pyquante2.graphics",
        "pyquante2.grid",
        "pyquante2.ints",
        "pyquante2.pt",
        "pyquante2.scf",
        "pyquante2.viewer",
    ],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
)
