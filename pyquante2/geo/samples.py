"""
A collection of molecules for testing and fun.

>>> h
[(1, 0.0, 0.0, 0.0)]

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 2.0 and later is covered by the GPL
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from pyquante2.geo.molecule import molecule

h = molecule([(1,0,0,0)],name='Hydrogen')

h2 = molecule([(1,  0.00000000,     0.00000000,     0.36628549),
               (1,  0.00000000,     0.00000000,    -0.36628549)],
              units='Angstrom',
              name='Hydrogen')

h2o = molecule([(8,   0.00000000,     0.00000000,     0.04851804),
                (1,   0.75300223,     0.00000000,    -0.51923377),
                (1,  -0.75300223,     0.00000000,    -0.51923377)],
               units='Angstrom',
               name='Water')

oh = molecule([(8,  0.00000000,     0.00000000,    -0.08687037),
               (1,  0.00000000,     0.00000000,     0.86464814)],
              units='Angstrom',
              multiplicity=2,
              name='Hydroxide')
he = molecule(atomlist = [(2,0,0,0)],name='Helium')
li = molecule(atomlist = [(3,0,0,0)], multiplicity=2,name='Lithium')
li_p = molecule(atomlist = [(3,0,0,0)],charge=1,name="Li+")
li_m = molecule(atomlist = [(3,0,0,0)],charge=-1,name="Li-")

lih = molecule([(3,    0.00000000,     0.00000000,    -0.53999756),
                (1,    0.00000000,     0.00000000,     1.08999756)],
               units='Angstrom',
               name="LiH")
             
co = molecule([(6,  0.00000000,     0.00000000,    -0.63546711),
               (8,  0.00000000,     0.00000000,     0.47832425)],
              units='Angstrom',
              name="CO")

ch4 = molecule([(6,   0.00000000,     0.00000000,     0.00000000),
                (1,   0.62558332,    -0.62558332,     0.62558332),
                (1,  -0.62558332,     0.62558332,     0.62558332),
                (1,   0.62558332,     0.62558332,    -0.62558332),
                (1,  -0.62558332,    -0.62558332,    -0.62558332)],
               units='Angstrom',
               name="CH4")

c6h6 = molecule([ (6,  0.98735329,     0.98735329,     0.00000000),
                  (6,  1.34874967,    -0.36139639,     0.00000000),
                  (6,  0.36139639,    -1.34874967,     0.00000000),
                  (6, -0.98735329,    -0.98735329,     0.00000000),
                  (6, -1.34874967,     0.36139639,     0.00000000),
                  (6, -0.36139639,     1.34874967,     0.00000000),
                  (1,  1.75551741,     1.75551741,     0.00000000),
                  (1,  2.39808138,    -0.64256397,     0.00000000),
                  (1,  0.64256397,    -2.39808138,     0.00000000),
                  (1, -1.75551741,    -1.75551741,     0.00000000),
                  (1, -2.39808138,     0.64256397,     0.00000000),
                  (1, -0.64256397,     2.39808138,     0.00000000)],
                units='Angstrom',
                name="Benzene")

if __name__ == '__main__':
    import doctest
    doctest.testmod()

