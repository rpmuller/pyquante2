#!/usr/bin/env python

from setuptools import setup
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
except ImportError:
    file_ext = "c"
    use_cython = False
else:
    import numpy as np

    file_ext = "pyx"
    use_cython = True


ext_modules = [
    Extension("pyquante2.cints.one", [f"cython/cone_wrap.{file_ext}", "cython/cints.c"]),
    Extension("pyquante2.cints.two", [f"cython/ctwo_wrap.{file_ext}", "cython/cints.c"]),
    Extension(
        "pyquante2.cints.hgp",
        [f"cython/chgp_wrap.{file_ext}", "cython/chgp.c", "cython/cints.c"],
    ),
    Extension("pyquante2.cints.rys", [f"cython/crys_wrap.{file_ext}", "cython/crys.c"]),
]

if use_cython:
    ext_modules.append(
        Extension(
            "pyquante2.cbecke", ["cython/cbecke.pyx"], include_dirs=[np.get_include()]
        )
    )
    ext_modules = cythonize(ext_modules)

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
    ext_modules=ext_modules,
)
