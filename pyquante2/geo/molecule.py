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
import numpy as np
from pyquante2 import settings
from pyquante2.geo.atom import atom
from pyquante2.utils import pairs,norm2

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

    def nuclear_repulsion(self):
        enuke = 0
        for i,ati in enumerate(self.atoms):
            for atj in self.atoms[:i]:
                enuke += ati.atno*atj.atno/np.sqrt(norm2(ati.r-atj.r))
        return enuke

    def nel(self):
        "Number of electrons of the molecule"
        return sum(atom.atno for atom in self) - self.charge
    
    def nocc(self,assert_closed=False):
        "Number of occupied orbitals"
        c,o = divmod(self.nel(),2)
        if assert_closed:
            assert o==0, "Molecule should be closed shell"
        return c+o
        


if __name__ == '__main__':
    import doctest
    doctest.testmod()
