"""
Create a molecule for use in pyquante

>>> h = molecule([(1,0,0,0)])
>>> h
Stoichiometry = H, Charge = 0, Multiplicity = 2
1 H     0.000000     0.000000     0.000000

>>> he = molecule([(2,0,0,0)])
>>> he
Stoichiometry = He, Charge = 0, Multiplicity = 1
2 He     0.000000     0.000000     0.000000
>>> he.nclosed(),he.nopen(),he.nup(),he.ndown()
(1, 0, 1, 1)

>>> he_trip = molecule([(2,0,0,0)], multiplicity=3)
>>> he_trip
Stoichiometry = He, Charge = 0, Multiplicity = 3
2 He     0.000000     0.000000     0.000000
>>> he_trip.nclosed(),he_trip.nopen(),he_trip.nup(),he_trip.ndown()
(0, 2, 2, 0)

"""
import numpy as np
from pyquante2 import settings
from pyquante2.geo.atom import atom
from pyquante2.utils import upairs

class molecule(object):
    """
    >>> from pyquante2.geo.samples import h2
    >>> np.isclose(h2.nuclear_repulsion(),0.722356)
    True
    >>> h2.nel()
    2
    >>> h2.nocc()
    1
    """
    def __init__(self,atomlist=None,**kwargs):
        self.atoms = []
        self.charge = int(kwargs.get('charge',settings.molecular_charge))
        self.multiplicity = kwargs.get('multiplicity')
        self.name = kwargs.get('name','pyquante2 molecule')

        self.units = kwargs.get('units',settings.units).lower()

        if atomlist is not None:
            for atuple in atomlist:
                self.atoms.append(atom(*atuple,units=self.units))
        if self.multiplicity is None:
            self.set_multiplicity()
        return

    def __getitem__(self,i): return self.atoms.__getitem__(i)
    def __len__(self): return self.atoms.__len__()

    def __repr__(self):
        lines = ["Stoichiometry = %s, Charge = %d, Multiplicity = %d" %\
                 (self.stoich(),self.charge,self.multiplicity)]
        lines.extend(repr(atom) for atom in self.atoms)
        return "\n".join(lines)

    def set_multiplicity(self):
        if self.nel() % 2:
            self.multiplicity = 2
        else:
            self.multiplicity = 1
        return

    def html(self,tablehead=True):
        import xml.etree.ElementTree as ET
        top = ET.Element("p")
        h2 = ET.SubElement(top,"h2")
        h2.text = self.name
        p = ET.SubElement(top,"p")
        p.text = "Stoichiometry = %s, Charge = %d, Multiplicity = %d" %\
                 (self.stoich(),self.charge,self.multiplicity)
        table = ET.SubElement(top,"table")
        if tablehead:
            tr = ET.SubElement(table,"tr")
            for item in ["#","Atno","Symbol","x","y","z"]:
                th = ET.SubElement(tr,"th")
                th.text = item
        for i,atom in enumerate(self.atoms):
            atom.html_row(table,i)
        return top

    def _repr_html_(self,tablehead=True):
        import xml.etree.ElementTree as ET
        top = ET.Element("html")
        top.append(self.html(tablehead=tablehead))
        return ET.tostring(top)

    def nuclear_repulsion(self):
        return sum(ati.atno*atj.atno/ati.distance(atj) for ati,atj in upairs(self))

    def nel(self):
        "Number of electrons of the molecule"
        return sum(atom.atno for atom in self) - self.charge
    
    def nocc(self): return sum(divmod(self.nel(),2))
    def nclosed(self):
        nc,ierr = divmod(self.nel()-self.nopen(),2)
        assert ierr == 0, "Error in molecule, multiplicity, nclosed, nopen"
        return nc
    def nopen(self): return self.multiplicity-1
    def nup(self): return self.nclosed()+self.nopen()
    def ndown(self): return self.nclosed()

    def xyz(self,title=None,fobj=None):
        """
        Output molecule in [xyz format](http://en.wikipedia.org/wiki/XYZ_file_format).
        """
        if title is None:
            title = self.name
        lines = ["%d" % len(self.atoms),"%s" % title]
        for atom in self.atoms:
            lines.append(atom.xyz())
        record = "\n".join(lines)
        if fobj:
            fobj.write(record)
        else:
            print(record)
        return

    def pyquante1(self,name="pyq2 molecule"):
        """
        Make a PyQuante1 Molecule object that can be passed into that program for
        testing/debugging purposes.
        """
        from PyQuante import Molecule
        atuples = [(a.atno,tuple(a.r)) for a in self.atoms]
        return Molecule(name,atuples,charge=self.charge,multiplicity=self.multiplicity)

    def bbox(self,padding=5.,BIG=1e12):
        xmin = ymin = zmin = BIG
        xmax = ymax = zmax = -BIG
        for atom in self.atoms:
            x,y,z = atom.r
            xmin = min(x,xmin)
            ymin = min(y,ymin)
            zmin = min(z,zmin)
            xmax = max(x,xmax)
            ymax = max(y,ymax)
            zmax = max(z,zmax)
        xmin,ymin,zmin = xmin-padding,ymin-padding,zmin-padding
        xmax,ymax,zmax = xmax+padding,ymax+padding,zmax+padding
        return xmin,xmax,ymin,ymax,zmin,zmax
        

    def stoich(self):
        """
        Generate a stoichiometry string for the molecule:
        >>> from pyquante2 import h2,h2o,c6h6
        >>> h2.stoich()
        'H2'
        >>> c6h6.stoich()
        'H6C6'
        """
        from collections import Counter
        from pyquante2.geo.elements import symbol
        cnt = Counter()
        for atom in self.atoms:
            cnt[atom.atno] += 1
        keys = sorted(cnt.keys())
        s = []
        for key in keys:
            if cnt[key] == 1:
                s.append(symbol[key])
            else:
                s.append("%s%d" % (symbol[key],cnt[key]))
        return "".join(s)

    def mass(self): return sum(at.mass() for at in self.atoms)
    def com(self): return sum(at.mass()*at.r for at in self.atoms)/self.mass()

    def center(self):
        r_com = self.com()
        for at in self.atoms:
            at.r = at.r - r_com
        return

    def bonds(self,scale=1.1):
        from pyquante2.utils import pairs
        bonds = []
        for i,ati in enumerate(self.atoms):
            for j in range(i):
                atj = self.atoms[j]
                if ati.distance(atj) < scale*(ati.radius()+atj.radius()):
                    bonds.append((i,j))
        return bonds

    def get_index(self,center,tol=1e-2):
        """\
        >>> he = molecule([(2,0,0,0)])
        >>> he.get_index((0,0,0))
        0
        """
        center = np.asarray(center,'d')
        for i,atom in enumerate(self.atoms):
            if np.linalg.norm(center-atom.r) < tol:
                return i
        return None

    def get_tag(self,center,tol=1e-2):
        """\
        >>> he = molecule([(2,0,0,0)])
        >>> he.get_tag((0,0,0))
        'He0'
        """
        index = self.get_index(center,tol)
        return self.atoms[index].tag(index)

def read_xyz(fname):
    with open(fname) as f:
        line = f.readline()
        nat = int(line.strip())
        comment = f.readline()
        lines = [f.readline() for i in range(nat)]
    return read_xyz_lines(lines)

def read_xyz_lines(lines,**kwargs):
    from pyquante2.geo.elements import sym2no
    from pyquante2.utils import parseline
    atuples = []
    for line in lines:
        sym,x,y,z = parseline(line,'sfff')
        atno = sym2no[sym]
        atuples.append((atno,x,y,z))
    return molecule(atuples,**kwargs)
                 
if __name__ == '__main__':
    import doctest
    doctest.testmod()
