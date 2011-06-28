"""
A collection of molecules for testing and fun.

>>> h
[(1, 0.0, 0.0, 0.0)]

Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from pyquante2.geo.molecule import molecule

h = molecule([(1,0,0,0)])

if __name__ == '__main__':
    import doctest
    doctest.testmod()

