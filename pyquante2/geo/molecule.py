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
from pyquante2.utils import upairs,norm2

class molecule:
    """
    >>> from pyquante2.geo.samples import h2
    >>> round(h2.nuclear_repulsion(), 6)
    0.72236
    >>> h2.nel()
    2
    >>> h2.nocc()
    1
    """
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
        return sum(ati.atno*atj.atno/np.sqrt(norm2(ati.r-atj.r)) for ati,atj in upairs(self))

    def nel(self):
        "Number of electrons of the molecule"
        return sum(atom.atno for atom in self) - self.charge
    
    def nocc(self): return sum(divmod(self.nel(),2))
    def nclosed(self): return self.nel()//2
    def nopen(self): return divmod(self.nel(),2)[1]
    def nup(self): return self.nocc()
    def ndown(self): return self.nclosed()
        
        


if __name__ == '__main__':
    import doctest
    doctest.testmod()
