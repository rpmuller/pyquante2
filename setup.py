from setuptools import setup, Extension
import numpy as np
import os

# Get the path to the directory containing your C headers
# This ensures the path is absolute and correct regardless of where the build starts
C_ROOT = os.path.abspath("src/pyquante2/cints/")

setup(
    ext_modules=[
        Extension("pyquante2.cints.one", 
                  ["src/pyquante2/cints/cone_wrap.c","src/pyquante2/cints/cints.c"],
                  include_dirs=[np.get_include(), C_ROOT]),
        Extension("pyquante2.cints.two", 
                  ["src/pyquante2/cints/ctwo_wrap.c", "src/pyquante2/cints/cints.c"],
                  include_dirs=[C_ROOT]),
        Extension("pyquante2.cints.hgp",
                  ["src/pyquante2/cints/chgp_wrap.c", "src/pyquante2/cints/chgp.c", "src/pyquante2/cints/cints.c"],
                  include_dirs=[C_ROOT]),
        Extension("pyquante2.cints.rys",
                  ["src/pyquante2/cints/crys_wrap.c", "src/pyquante2/cints/crys.c"]),
        Extension("pyquante2.cbecke",
                  ["src/pyquante2/cints/cbecke.c"], 
                  include_dirs=[np.get_include()]),
        ],
)
