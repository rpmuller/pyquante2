"""
Create a molecule for use in pyquante
>>> h = molecule([(1,0,0,0)])
>>> h
[(1, 0.0, 0.0, 0.0)]

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 2.0 and later is covered by the GPL
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from pyquante2 import settings
from pyquante2.geo.atom import atom

class molecule:
    def __init__(self,atomlist=[],**kwargs):
        self.atoms = []
        self.charge = int(kwargs.get('charge',settings.molecular_charge))
        self.multiplicity = int(kwargs.get('multiplicity',settings.spin_multiplicity))

        self.units = kwargs.get('units',settings.units).lower()

        if atomlist:
            for atuple in atomlist:
                self.atoms.append(atom(*atuple,units=self.units))
        return

    def __repr__(self): return repr(self.atoms)

    def __getitem__(self,i): return self.atoms.__getitem__(i)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
