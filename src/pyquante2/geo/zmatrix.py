"""\
Utilities to read/write zmatrices.

Todo:
- Write a function to create a zmatrix from a molecule, possibly given
  a set of internal coordinates, possibly not.
"""

from math import pi,cos,sin
from numpy import array,cross,dot
from numpy.linalg import norm

def parse_zmatrix(multistring):
    """\
    Parse a multiline string into zmatrix format.
    >>> parse_zmatrix("H")
    [['H']]

    >>> parse_zmatrix('H\\nH 1 0.7')
    [['H'], ['H', 1, 0.7]]
    """
    parsers = [str,int,float,int,float,int,float]
    lines = []
    for line in multistring.splitlines():
        words = line.split()
        tokens = []
        lines.append(tokens)
        for i,word in enumerate(words):
            tokens.append(parsers[i](word))
    return lines

def unpack_zmat_line(words,index=0):
    def radians(ang): return ang*pi/180.0
    
    sym, I,R, J,theta, K,phi = 0, 0,0, 0,0, 0,0
    sym = words[0]
    if index > 0:
        I = words[1]
        R = words[2]
    if index > 1:
        J = words[3]
        theta = radians(words[4])
    if index > 2:
        K = words[5]
        phi = radians(words[6])
    return sym,I,R,J,theta,K,phi

def z2xyz(geo):
    """
    Convert geometry from zmatrix coordinates to xyz coordinates.
    
    geo is a list of tuples containing the zmatrix information.
    >>> z2xyz([['H'], ['H', 1, 0.7]])
    [['H', 0, 0, 0], ['H', 0.7, 0, 0]]

    J. Parsons, et al., J Comput Chem 26, 1063 (2005) is a good ref for this.
    """
    xyzs = []
    for i,atom in enumerate(geo):
        sym,I,R,J,theta,K,phi = unpack_zmat_line(atom,i)

        if i == 0:
            xyzs.append([sym,0,0,0])
        elif i == 1:
            xyzs.append([sym,R,0,0])
        elif i == 2:
            B = array(xyzs[I-1][1:])
            BA = array([-cos(theta),sin(theta),0])
            xyz = B+R*BA
            xyzs.append([sym]+xyz.tolist())
        else:
            B = array(xyzs[I-1][1:])
            C = array(xyzs[J-1][1:])
            D = array(xyzs[K-1][1:])
            CD = D-C
            BC = C-B
            n1 = cross(CD,BC)
            n1 /= norm(n1)
            n2 = rotate(n1,BC,phi)
            n2 /= norm(n2)
            BA = rotate(BC,n2,-theta)
            BA /= norm(BA)
            xyz = B+R*BA
            xyzs.append([sym]+xyz.tolist())
    return xyzs

def zmatrix_tomolecule(zmt):
    "Create a pyquante2 molecule from a zmatrix"
    from pyquante2.geo.elements import sym2no
    from pyquante2.geo.molecule import molecule
    atuples = [(sym2atno[sym],x,y,z) for sym,x,y,z in z2xyz(geo)]
    return molecule(atuples)

def zmatrix_tostring(zmt):
    lines = []
    for i,z in enumerate(zmt):
        if i == 0:
            lines.append("%4s" % z[0])
        elif i == 1:
            lines.append("%4s %4i %10.4f" % tuple(z))
        elif i == 2:
            lines.append("%4s %4i %10.4f %4i %10.4f" % tuple(z))
        else:
            lines.append("%4s %4i %10.4f %4i %10.4f %4i %10.4f" % tuple(z))
    return "\n".join(lines)

def rotate(v,k,theta):
    """[Euler-Rogrigues rotation formula](http://en.wikipedia.org/wiki/Rodrigues'_rotation_formula)
    """
    c = cos(theta)
    s = sin(theta)
    return v*c + cross(k,v)*s + k*dot(k,v)*(1-c)

def isclose(x1,x2,tol=1e-8): return abs(x1-x2) < tol

def cartesians_equal(c1,c2,tol=1e-8,verbose=True):
    for at1,at2 in zip(c1,c2):
        if not (at1[0] == at2[0] and\
                isclose(at1[1],at2[1],tol) and\
                isclose(at1[2],at2[2],tol) and\
                isclose(at1[3],at2[3],tol)):
            return False
    return True
           
def test():
    assert parse_zmatrix("H") == [['H']]
    assert parse_zmatrix("H\nH 1 0.7") == [['H'], ['H', 1, 0.7]]
    assert parse_zmatrix("O\nH 1 1.0\nH 1 1.0 2 90.0") == [['O'],['H',1,1.0],['H',1,1.0,2,90.0]]
    
    assert cartesians_equal(z2xyz([['H']]),
                            [['H',0.0,0.0,0.0]])
    assert cartesians_equal(z2xyz([['H'], ['H', 1, 0.7]]),
                            [['H',0.0,0.0,0.0],['H',0.7,0.0,0.0]])
    assert cartesians_equal(z2xyz([['O'],['H',1,1.0],['H',1,1.0,2,90.0]]),
                            [['O',0.0,0.0,0.0],['H',1.0,0.0,0.0],['H',0.0,1.0,0.0]])
    assert cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,180.]]),
                            [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',2,1,0]])
    assert cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,0.]]),
                            [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',0,1,0]])
    assert cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,90.]]),
                            [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',1,1,-1]])
    assert cartesians_equal(z2xyz([['H'],['O',1,1],['O',2,1,1,90.],['H',3,1,2,90.,1,-90.]]),
                            [['H',0,0,0],['O',1,0,0],['O',1,1,0],['H',1,1,1]])
    return

def interp(x1,x2,amt): return amt*x1+(1-amt)*x2

def simple_zmatrix_interp(z1,z2,amt):
    zout = []
    for i,(z1,z2) in enumerate(zip(z1,z2)):
        if i == 0:
            s1 = z1[0]
            s2 = z2[0]
            assert s1==s2
            zout.append([s1])
        elif i == 1:
            s1,i1,r1 = z1
            s2,i2,r2 = z2
            assert s1==s2
            assert i1==i2
            zout.append([s1,i1,interp(r1,r2,amt)])
        elif i == 2:
            s1,i1,r1,j1,th1 = z1
            s2,i2,r2,j2,th2 = z2
            assert s1==s2
            assert i1==i2
            assert j1==j2
            zout.append([s1,i1,interp(r1,r2,amt),j1,interp(th1,th2,amt)])
        else:
            s1,i1,r1,j1,th1,k1,phi1 = z1
            s2,i2,r2,j2,th2,k2,phi2 = z2
            assert s1==s2
            assert i1==i2
            assert j1==j2
            assert k1==k2
            zout.append([s1,i1,interp(r1,r2,amt),j1,interp(th1,th2,amt),k1,interp(phi1,phi2,amt)])
    return zout
