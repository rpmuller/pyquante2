"""
CML module

This module contains some general functions I use for parsing CML
files. It currently uses an intermediate representation to store the data
rather than pyquante2 data structures.
"""

import io
import numpy
import xml.etree.ElementTree as ET
from pyquante2.geo.elements import sym2no,mass,rvdw

# For testing purposes
waters_cml = io.StringIO("""\
<cml>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.0" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.1" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.2" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.3" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.4" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.5" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.6" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.7" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.8" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
<molecule>
  <name>Water</name>
  <atomArray>
    <atom id="a1" x3="0.0" y3="0.0" z3="0.0" elementType="O"/>
    <atom id="a2" x3="1.0" y3="0.0" z3="0.0" elementType="H"/>
    <atom id="a3" x3="0.0" y3="1.9" z3="0.0" elementType="H"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a1 a3" order="1"/>
  </bondArray>
</molecule>
</cml>""")

def add_atomarray(mol): return subelement(mol,"atomArray")

def add_atomarray_and_atoms(mol,atoms,units):
    atomarray = add_atomarray(mol)
    for i,(sym,x,y,z) in enumerate(atoms):
        atom = add_atom(atomarray,sym,x,y,z,units,i)
    return atomarray

def add_bond(bondarray,idi,idj):
    bond = subelement(bondarray,'bond',atomRefs2="%s %s" % (idi,idj),
                      order="1")
    return bond

def add_bondarray(mol):
    bondarray = subelement(mol,'bondArray')
    return bondarray

def add_bondarray_and_bonds(mol,bonds):
    bondarray = add_bondarray(mol)
    for idi,idj in bonds:
        bond = add_bond(bondarray,idi,idj)
    return bondarray

def add_cartatom(root,sym,x,y,z,units,i):
    atomel = subelement(root,'atom',id="a%d"%(i+1),elementType=sym,
                        x3=str(x),y3=str(y),z3=str(z))
    return atomel

def add_mol(parent,atoms,uc=None,units='cartesian',**kwargs):
    mol = subelement(parent,'molecule')
    if uc is not None:
        add_crystal(mol,uc)
    atomarray = add_atomarray_and_atoms(mol,atoms,units)
    bonds = find_bonds_in_mol(mol)
    if bonds:
        bondarray = add_bondarray_and_bonds(mol,bonds)
    return mol

def add_mol(parent,atoms,uc=None,units='cartesian',**kwargs):
    mol = subelement(parent,'molecule')
    if uc is not None:
        add_crystal(mol,uc)
    atomarray = add_atomarray_and_atoms(mol,atoms,units)
    bonds = find_bonds_in_mol(mol)
    if bonds:
        bondarray = add_bondarray_and_bonds(mol,bonds)
    return mol

def atomtuples(mol):
    "Convert a molecule subtree to sym,x,y,z atomtuple"
    return [atomtuple(atom) for atom in get_atoms(mol)]
        
def atomtuple(atom):
    "Convert an atom subtree to sym,x,y,z atomtuple"
    x,y,z = get_xyz(atom)
    return get_sym(atom),x,y,z

def cart_distance(dx,dy,dz):
    return math.sqrt(dx*dx+dy*dy+dz*dz)

def center_of_mass(molecule):
    """
    Given a molecule rooted at *molecule*, compute the center of mass.
    For this routine we're going to assume we have non-fractional coordinates.
    """
    xcom=ycom=zcom=0
    totm = 0
    for atom in get_atoms(molecule):
        m = get_mass(atom)
        x,y,z = get_xyz(atom)
        xcom += m*x
        ycom += m*y
        zcom += m*z
        totm += m
    xcom /= totm
    ycom /= totm
    zcom /= totm
    return xcom,ycom,zcom

def find_bonds_in_mol(mol,bondscale=1.1):
    bonds = []
    atoms = get_atoms(mol)
    nat = len(atoms)
    for i in range(nat):
        ati = atoms[i]
        ri0 = get_vdw_radius(ati)
        xi,yi,zi = get_xyz(ati)
        idi = get_id(ati)
        for j in range(i):
            atj = atoms[j]
            rj0 = get_vdw_radius(atj)
            xj,yj,zj = get_xyz(atj)
            idj = get_id(atj)
                
            dx,dy,dz = xi-xj,yi-yj,zi-zj
            r0 = (ri0+rj0)/2
            r = cart_distance(dx,dy,dz)
            if r < bondscale*r0:
                bonds.append((idi,idj))
    return bonds

def getattrib(tag,atname,default=None): return tag.attrib.get(atname,default)

def getattribf(tag,atname,default=None):
    return float(getattrib(tag,atname,default))

def gettext(tag): return tag.text

def gettextf(tag): return float(gettext(tag))

def get_atoms(root):
    return get_records_with_tag(root,'atom')

def get_atno(atom): return sym2no[get_sym(atom)]

def get_id(atom):
    return getattrib(atom,"id")

def get_mass(atom): return mass[get_atno(atom)]

def get_molecules(root):
    """Return a list of trees rooted at all of the molecules in
    the tree rooted at root"""
    mols = get_records_with_tag(root,'molecule')
    if not mols: mols = [root] # Assume tree consists of a molecule at the root
    return mols

def get_records_with_tag(inp,tag):
    """
    Given a filename or an element tree in *inp*, return all records
    matching *tag*
    """
    if isinstance(inp, str):  # Assume inp is a filename
        inp = read_cml(inp)
    return inp.findall(".//%s" % tag)

def get_sym(atom):
    return getattrib(atom,"elementType")

def get_vdw_radius(atom): return rvdw[get_atno(atom)]

def get_xyz(atom):
    return get_xyzc(atom)

def get_xyzc(atom):
    return getattribf(atom,"x3"),getattribf(atom,"y3"),\
           getattribf(atom,"z3")

def hasattrib(tag,atname): return tag in atname.attrib

def natoms(mol):
    return len(get_atoms(mol))

def randomindex(lim=10**8):
    import random
    return random.randrange(lim)
    
def read_cml(fname):
    return ET.parse(fname)

def read(fname):
    tree = ET.parse(fname)
    #ET.dump(tree)
    return tree

def remove_bonds(root):
    for mol in get_molecules(root):
        for ba in get_records_with_tag(mol,'bondArray'):
            mol.remove(ba)
    return

def setattrib(tag,atname,val): tag.attrib[atname]=str(val)

def stoichiometry(mol):
    stdict = {}
    for atom in get_atoms(mol):
        sym = get_sym(atom)
        if sym in stdict:
            stdict[sym] += 1
        else:
            stdict[sym] = 1
    return stoichiometry_formatter(stdict)

def stoichiometry_formatter(stdict):
    l = []
    for sym in stdict:
        l.append("%s%d" % (sym,stdict[sym]))
    return "".join(l)

def subelement(parent,name,text=None,**kwargs):
    # parent can be none!
    el = ET.Element(name)
    if parent is not None: parent.append(el)
    if text:
        el.text = str(text)
    for key in kwargs:
        el.set(key,str(kwargs[key]))
    return el

def test():
    waters = read(waters_cml)
    for mol in get_molecules(waters):
        print("Mol: ")
        for atom in get_atoms(mol):
            ET.dump(atom)
        print(center_of_mass(mol))
    return

def tocml(atoms,uc=None,units='cartesian',**kwargs):
    mol = add_mol(None,atoms,uc,units,**kwargs)
    tree = ET.ElementTree(mol)
    return tree

def tocml_multigeo(geos,uc=None,units='cartesian',**kwargs):
    cml = ET.Element('cml')
    for geo in geos:
        if len(geo):
            mol = add_mol(cml,geo,uc,units,**kwargs)
    tree = ET.ElementTree(cml)
    return tree

def translate(molecule,dx,dy,dz):
    """
    Given a molecule rooted at *molecule*, translate all atoms by
    dx,dy,dz
    """
    for atom in get_atoms(molecule):
        translate_atom(atom,dx,dy,dz)
    return

def translate_atom(atom,dx,dy,dz):
    setattrib(atom,"x3",getattribf(atom,"x3")+dx)
    setattrib(atom,"y3",getattribf(atom,"y3")+dy)
    setattrib(atom,"z3",getattribf(atom,"z3")+dz)
    return

def translate_com(tree):
    xcom,ycom,zcom = center_of_mass(tree)
    translate(tree,-xcom,-ycom,-zcom)
    return

def write_cml(fname,atoms,uc=None,units='cartesian',**kwargs):
    tree = tocml(atoms,uc,units,**kwargs)
    tree.write(fname)
    if kwargs.get('verbose'):
        ET.dump(tree)
    return

def write_cml_multigeo(fname,geos,uc=None,units='cartesian',**kwargs):
    tree = tocml_multigeo(geos,uc,units,**kwargs)
    tree.write(fname)
    if kwargs.get('verbose'): ET.dump(tree)
    return

if __name__ == '__main__':
    test()

