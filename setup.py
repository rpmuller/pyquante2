from setuptools import setup, Extension
import numpy as np

setup(
    ext_modules=[
        Extension("pyquante2.cints.one", 
                  ["src/pyquante2/cints/cone_wrap.c","src/pyquante2/cints/cints.c"]),
        Extension("pyquante2.cints.two", 
                  ["src/pyquante2/cints/ctwo_wrap.c", "src/pyquante2/cints/cints.c"]),
        Extension("pyquante2.cints.hgp",
                  ["src/pyquante2/cints/chgp_wrap.c", "src/pyquante2/cints/chgp.c", "src/pyquante2/cints/cints.c"]),
        Extension("pyquante2.cints.rys",
                  ["src/pyquante2/cints/crys_wrap.c", "src/pyquante2/cints/crys.c"]),
        Extension("pyquante2.cbecke",
                  ["src/pyquante2/cints/cbecke.c"], 
                  include_dirs=[np.get_include()]),
        ],
)
