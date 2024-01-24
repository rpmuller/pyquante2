from pathlib import Path
from setuptools import setup
from setuptools.extension import Extension

import numpy as np

try:
    from Cython.Build import cythonize
except ImportError:
    FILE_EXT = "c"
    USE_CYTHON = False
else:
    FILE_EXT = "pyx"
    USE_CYTHON = True


_NUMPY_MACROS = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
ext_modules = [
    Extension(
        "pyquante2.cints.one",
        [f"cython/cone_wrap.{FILE_EXT}", "cython/cints.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.two",
        [f"cython/ctwo_wrap.{FILE_EXT}", "cython/cints.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.hgp",
        [f"cython/chgp_wrap.{FILE_EXT}", "cython/chgp.c", "cython/cints.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.rys",
        [f"cython/crys_wrap.{FILE_EXT}", "cython/crys.c"],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cbecke",
        [f"cython/cbecke.{FILE_EXT}"],
        define_macros=_NUMPY_MACROS,
        include_dirs=[np.get_include()],
    ),
]

if USE_CYTHON:
    ext_modules = cythonize(ext_modules, annotate=True)

# Make intermediate directories so that `setup.py build_ext -i` works.
_PWD = Path(".").resolve()
for ext_module in ext_modules:
    # The last component is the compiled library itself.
    path = _PWD.joinpath(*ext_module.name.split(".")[:-1])
    path.mkdir(parents=True, exist_ok=True)

setup(
    name="pyquante2",
    ext_modules=ext_modules,
)
