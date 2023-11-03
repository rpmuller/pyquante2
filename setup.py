from pathlib import Path
from setuptools import setup
from setuptools.extension import Extension

try:
    from Cython.Build import cythonize
except ImportError:
    FILE_EXT = "c"
    USE_CYTHON = False
else:
    import numpy as np

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
]

if USE_CYTHON:
    ext_modules.append(
        Extension(
            "pyquante2.cbecke",
            ["cython/cbecke.pyx"],
            define_macros=_NUMPY_MACROS,
            include_dirs=[np.get_include()],
        )
    )
    ext_modules = cythonize(ext_modules, annotate=True)

# Make intermediate directories so that `setup.py build_ext -i` works.
_PWD = Path(".").resolve()
for ext_module in ext_modules:
    # The last component is the compiled library itself.
    path = _PWD.joinpath(*ext_module.name.split(".")[:-1])
    path.mkdir(parents=True, exist_ok=True)

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
