import logging
from pyquante2.grid.grid import grid
from pyquante2.ints.integrals import onee_integrals,twoe_integrals
from pyquante2.utils import trace2, geigh
from pyquante2.scf.iterators import SCFIterator,USCFIterator,AveragingIterator,ROHFIterator
import numpy as np

logger = logging.getLogger(__name__)

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
    >>> np.isclose(h2_rhf.energy,-1.11709942949)
    True

    >>> ens = h2_rhf.converge(AveragingIterator,maxiters=100)
    >>> np.isclose(h2_rhf.energy,-1.11709325545)
    True

    >>> ens = h2_rhf.converge(SCFIterator,maxiters=1)
    >>> np.isclose(h2_rhf.energy,-1.1170994294946217)
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
        self.energy += trace2(H,D)
        E,c = geigh(H,self.i1.S)
        self.orbe = E
        self.orbs = c
        return c

class dft(rhf):
    "Hamiltonian for DFT calculations. Adds a grid to RHF iterator."
    def __init__(self,geo,bfs,xcname='lda',verbose=False):
        rhf.__init__(self,geo,bfs)
        self.grid = grid(geo)
        self.xcname = xcname
        self.grid.setbfamps(bfs)
        self.verbose = verbose
        return

    def update(self,D):
        from pyquante2.dft.dft import get_xc
        E0 = self.geo.nuclear_repulsion()
        self.energy = E0
        H = self.i1.T + self.i1.V
        E1 = 2*trace2(H,D)
        self.energy += E1

        J = self.i2.get_j(D)
        Ej = 2*trace2(J,D)

        # The 0.5 before the D comes from making the alpha density from the total density
        Exc,Vxc = get_xc(self.grid,0.5*D,xcname=self.xcname)

        H = H + 2*J + Vxc
        self.energy += Ej+Exc
        
        if self.verbose: logger.info(f"Energy: {self.energy}, E1: {E1}, Ej: {Ej}, Exc: {Exc}, E0: {E0}")
        E,c = geigh(H,self.i1.S)
        self.orbe = E
        self.orbs = c
        return c

class uhf(hamiltonian):
    """
    >>> from pyquante2.geo.samples import li
    >>> from pyquante2.basis.basisset import basisset
    >>> from pyquante2.scf.iterators import USCFIterator
    >>> bfs = basisset(li,'6-31G**')
    >>> solver = uhf(li,bfs)
    >>> ens = solver.converge(USCFIterator)
    >>> np.isclose(solver.energy,-7.4313707537)
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

class rohf(hamiltonian):
    """Hamiltonian for ROHF calculations. This is the simple version from ???,
    rather than WAG's version that also does GVB.

    >>> from pyquante2.geo.samples import he, he_triplet
    >>> from pyquante2.basis.basisset import basisset
    >>> bfs = basisset(he,'6-31G**')
    >>> he1 = rohf(he,bfs)
    >>> ens = he1.converge()
    >>> np.isclose(he1.energy,-2.855160702)
    True
    >>> he3 = rohf(he_triplet,bfs)
    >>> ens = he3.converge()
    >>> np.isclose(he3.energy,-1.3993077765340005)
    True
    """
    name = 'ROHF'
    
    def converge(self,iterator=ROHFIterator,**kwargs):
        return hamiltonian.converge(self,iterator,**kwargs)

    def update(self,Da,Db,orbs):
        from pyquante2.utils import ao2mo
        nalpha = self.geo.nup()
        nbeta = self.geo.ndown()
        norbs = len(orbs) # Da.shape[0]
        
        E0 = self.geo.nuclear_repulsion()
        h = self.i1.T + self.i1.V
        E1 = 0.5*trace2(Da+Db,h)
        Ja,Ka = self.i2.get_j(Da),self.i2.get_k(Da)
        Jb,Kb = self.i2.get_j(Db),self.i2.get_k(Db)
        Fa = h + Ja + Jb - Ka
        Fb = h + Ja + Jb - Kb
        E2 = 0.5*(trace2(Fa,Da)+trace2(Fb,Db))
        self.energy = E0+E1+E2
        #print (self.energy,E1,E2,E0)

        Fa = ao2mo(Fa,orbs)
        Fb = ao2mo(Fb,orbs)

        # Building the approximate Fock matrices in the MO basis
        F = 0.5*(Fa+Fb)
        K = Fb-Fa

        # The Fock matrix now looks like
        #      F-K    |  F + K/2  |    F
        #   ---------------------------------
        #    F + K/2  |     F     |  F - K/2
        #   ---------------------------------
        #       F     |  F - K/2  |  F + K

        # Make explicit slice objects to simplify this
        do = slice(0,nbeta)
        so = slice(nbeta,nalpha)
        uo = slice(nalpha,norbs)
        F[do,do] -= K[do,do]
        F[uo,uo] += K[uo,uo]
        F[do,so] += 0.5*K[do,so]
        F[so,do] += 0.5*K[so,do]
        F[so,uo] -= 0.5*K[so,uo]
        F[uo,so] -= 0.5*K[uo,so]

        E,cmo = np.linalg.eigh(F)
        c = np.dot(orbs,cmo)

        self.orbe = E
        self.orbs = c

        return c

        
if __name__ == '__main__':
    import doctest; doctest.testmod()
