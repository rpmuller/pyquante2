from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

_NUMPY_MACROS = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]

ext_modules = [
    Extension(
        "pyquante2.cints.one",
        ["cython/cone_wrap.pyx", "cython/cints.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.two",
        ["cython/ctwo_wrap.pyx", "cython/cints.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.hgp",
        ["cython/chgp_wrap.pyx", "cython/chgp.c", "cython/cints.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.rys",
        ["cython/crys_wrap.pyx", "cython/crys.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cbecke",
        ["cython/cbecke.pyx"],
        define_macros=_NUMPY_MACROS,
        include_dirs=[np.get_include()],
    ),
]

ext_modules = cythonize(ext_modules, annotate=True)

setup(
    name="pyquante2",
    ext_modules=ext_modules,
)

