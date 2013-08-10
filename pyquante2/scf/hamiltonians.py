from pyquante2.grid.grid import grid
from pyquante2.ints.integrals import onee_integrals,twoe_integrals
from pyquante2.utils import trace2, geigh, isnear
from pyquante2.scf.iterators import SCFIterator,USCFIterator,AveragingIterator
import numpy as np

class hamiltonian(object):
    name = 'abstract'
    def __init__(self,geo,bfs):
        self.geo = geo
        self.bfs = bfs
        self.i1 = onee_integrals(bfs,geo)
        self.i2 = twoe_integrals(bfs)
        self.energies = []
        self.energy = 0
        self.converged = False

    def __repr__(self):
        lines = ["%s Hamiltonian" % self.name]
        lines.append(str(self.geo))
        lines.append("Basis set: %s, Nbf: %d" %  (self.bfs.name,len(self.bfs)))
        lines.append("Status: Converged = %s" % self.converged)
        for i,E in enumerate(self.energies):
            lines.append("%d  %.5f" % (i,E))
        return "\n".join(lines)

    def _repr_html_(self):
        import xml.etree.ElementTree as ET
        top = ET.Element("html")
        h2 = ET.SubElement(top,"h2")
        h2.text = "%s Hamiltonian" % self.name
        top.append(self.geo.html())
        p = ET.SubElement(top,"p")
        p.text = "Basis set: %s, Nbf: %d" % (self.bfs.name,len(self.bfs))
        p = ET.SubElement(top,"p")
        p.text = "Status: Converged=%s" % self.converged
        if self.energies:
            table = ET.SubElement(top,"table")
            tr = ET.SubElement(table,"tr")
            for heading in ["#","Energy"]:
                td = ET.SubElement(tr,"th")
                td.text = heading
            for i,energy in enumerate(self.energies):
                tr = ET.SubElement(table,"tr")
                td = ET.SubElement(table,"td")
                td.text = str(i)
                td = ET.SubElement(table,"td")
                td.text = "%.5f" % energy
        return ET.tostring(top)

    def converge(self,iterator=SCFIterator,**kwargs):
        converger = iterator(self,**kwargs)
        self.energies = []
        for en in converger:
            self.energies.append(en)
        self.converged = converger.converged
        return self.energies

    def update(self,*args,**kwargs): raise Exception("Unimplemented")

class rhf(hamiltonian):
    """
    >>> from pyquante2.geo.samples import h2
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(h2,'sto3g')
    >>> h2_rhf = rhf(h2,bfs)
    >>> ens = h2_rhf.converge(SCFIterator)
    >>> isnear(h2_rhf.energy,-1.11709942949)
    True

    >>> ens = h2_rhf.converge(AveragingIterator,maxiters=100)
    >>> isnear(h2_rhf.energy,-1.11709325545)
    True

    >>> ens = h2_rhf.converge(SCFIterator,maxiters=1)
    >>> isnear(h2_rhf.energy,0.485554)
    True
    >>> h2_rhf.converged
    False
    """
    name = 'RHF'

    def update(self,D):
        self.energy = self.geo.nuclear_repulsion()
        H = self.i1.T + self.i1.V
        self.energy += trace2(H,D)

        JK = self.i2.get_2jk(D)
        H = H + JK
        #print(H)
        self.energy += trace2(H,D)
        E,c = geigh(H,self.i1.S)
        self.orbe = E
        self.orbs = c
        return c

class rdft(rhf):
    "Hamiltonian for DFT calculations. Adds a grid to RHF iterator."
    def __init__(self,geo,bfs):
        rhf.__init__(self,geo,bfs)
        self.grid = grid(geo)
        # make grid here.

    def update(self,D):
        self.energy = self.geo.nuclear_repulsion()
        H = self.i1.T + self.i1.V
        self.energy += trace2(H,D)

        J = self.i2.get_2j(D)
        H = H + J

        # XC = ???
        # H = H + XC

        self.energy += trace2(H,D)
        E,c = geigh(H,self.i1.S)
        self.orbe = E
        self.orbs = c
        return c

class rohf(rhf):
    """Hamiltonian for ROHF calculations. Adds shells information
    >>> from pyquante2.geo.samples import h2
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(h2,'sto3g')
    >>> h2_singlet = rohf(h2,bfs,[1],[1])
    >>> h2_triplet = rohf(h2,bfs,[1,1],[0.5,0.5])
    """
    def __init__(self,geo,bfs,norbsh=[],fi=[]):
        rhf.__init__(self,geo,bfs)
        self.norbsh = norbsh
        self.fi = fi

        
class uhf(hamiltonian):
    """
    >>> from pyquante2.geo.samples import oh
    >>> from pyquante2.basis.basisset import basisset
    >>> from pyquante2.scf.iterators import USCFIterator
    >>> bfs = basisset(oh,'sto3g')
    >>> solver = uhf(oh,bfs)
    >>> ens = solver.converge(USCFIterator)
    >>> isnear(solver.energy,-74.146669)
    True
    """
    name = 'UHF'

    def converge(self,iterator=USCFIterator,**kwargs):
        return hamiltonian.converge(self,iterator,**kwargs)

    def update(self,Da,Db):
        self.energy = self.geo.nuclear_repulsion()
        h = self.i1.T + self.i1.V
        self.energy += trace2(Da+Db,h)/2
        Ja,Ka = self.i2.get_j(Da),self.i2.get_k(Da)
        Jb,Kb = self.i2.get_j(Db),self.i2.get_k(Db)
        Fa = h + Ja + Jb - Ka
        Fb = h + Ja + Jb - Kb
        orbea,ca = geigh(Fa,self.i1.S)
        orbeb,cb = geigh(Fb,self.i1.S)
        self.energy += trace2(Fa,Da)/2 + trace2(Fb,Db)/2
        self.orbea = orbea
        self.orbsa = ca
        self.orbeb = orbeb
        self.orbsb = cb
        return ca,cb

if __name__ == '__main__':
    import doctest; doctest.testmod()
