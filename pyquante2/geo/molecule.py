"""
Create a molecule for use in pyquante
>>> h = molecule([(1,0,0,0)])
>>> h
[(1, 0.0, 0.0, 0.0)]

Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
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
        if atomlist:
            self.add_atuples(atomlist)
        return

    def __repr__(self): return repr(self.atoms)

    def __getitem__(self,i): return self.atoms.__getitem__(i)

    def add_atuples(self,atuples):
        for atuple in atuples:
            self.add_atuple(atuple)
        return

    def add_atuple(self,atuple): self.atoms.append(atom(*atuple))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
