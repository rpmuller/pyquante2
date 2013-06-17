from pyquante2.ints.integrals import onee_integrals,twoe_integrals
from pyquante2.utils import trace2,geigh
from pyquante2.scf.iterators import simple,usimple
import numpy as np

class hamiltonian:
    name = 'abstract'
    def __init__(self,geo,bfs):
        self.geo = geo
        self.bfs = bfs
        self.i1 = onee_integrals(bfs,geo)
        self.i2 = twoe_integrals(bfs)
        self.energies = []
        self.converged = False

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

class rhf(hamiltonian):
    """
    >>> from pyquante2.geo.samples import h2
    >>> from pyquante2.basis.basisset import basisset
    >>> from pyquante2.scf.iterators import simple,averaging
    >>> bfs = basisset(h2,'sto3g')
    >>> h2_rhf = rhf(h2,bfs)
    >>> ens = h2_rhf.converge(simple)
    >>> round(h2_rhf.energy,6)
    -1.1171
    """
    name = 'RHF'

    def converge(self,iterator=simple,**kwargs):
        self.energies = list(iterator(self,**kwargs))
        # Need something that checks for convergence, rather than just max iterations
        self.converged = True
        return self.energies

    def update(self,D):
        self.energy = self.geo.nuclear_repulsion()
        H = self.i1.T + self.i1.V
        self.energy += trace2(H,D)

        JK = self.i2.get_2jk(D)
        H = H + JK
        self.energy += trace2(H,D)
        E,c = geigh(H,self.i1.S)
        return c
        
class uhf(hamiltonian):
    name = 'UHF'

    def converge(self,iterator=usimple,**kwargs):
        self.energies = list(iterator(self,**kwargs))
        # Need something that checks for convergence, rather than just max iterations
        self.converged = True
        return self.energies

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
        return ca,cb

if __name__ == '__main__':
    import doctest; doctest.testmod()
