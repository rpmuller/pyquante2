import os.path
import sys
from os import getenv
from pathlib import Path
from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

import numpy as np

CYTHON_DISABLED = getenv("PYQUANTE_CYTHON_DISABLED", "False").lower() in (
    "true",
    "1",
    "t",
)

try:
    from Cython.Build import cythonize

    FILE_EXT = "pyx"
    CYTHON_AVAILABLE = True
    print("Cython available")
except ImportError:
    FILE_EXT = "c"
    CYTHON_AVAILABLE = False
    print("Cython not available")

if CYTHON_AVAILABLE and not CYTHON_DISABLED:
    print("Cython available and requested")
    USE_CYTHON = True
else:
    USE_CYTHON = False


_NUMPY_MACROS = [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
ext_modules = [
    Extension(
        "pyquante2.cints.one",
        [
            os.path.join("cython", f"cone_wrap.{FILE_EXT}"),
            os.path.join("cython", "cints.c"),
        ],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.two",
        [
            os.path.join("cython", f"ctwo_wrap.{FILE_EXT}"),
            os.path.join("cython", "cints.c"),
        ],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.hgp",
        [
            os.path.join("cython", f"chgp_wrap.{FILE_EXT}"),
            os.path.join("cython", "chgp.c"),
            os.path.join("cython", "cints.c"),
        ],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cints.rys",
        [
            os.path.join("cython", f"crys_wrap.{FILE_EXT}"),
            os.path.join("cython", "crys.c"),
        ],
        define_macros=_NUMPY_MACROS,
    ),
    Extension(
        "pyquante2.cbecke",
        [os.path.join("cython", f"cbecke.{FILE_EXT}")],
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


class optional_build_ext(build_ext):
    """build_ext that allows installation to continue if Extensions can't be built.

    This is useful for when there is code logic that will instead route to a
    pure Python implementation.
    """

    def build_extension(self, ext):
        self.couldnt_build = False
        try:
            super().build_extension(ext)
        except:
            self.couldnt_build = True
            print(
                "Couldn't compile extension, will fall back to pure Python implementations",
                file=sys.stderr,
            )

    def copy_extensions_to_source(self):
        if not self.couldnt_build:
            super().copy_extensions_to_source()


setup(
    name="pyquante2",
    ext_modules=ext_modules,
    cmdclass={"build_ext": optional_build_ext},
)
