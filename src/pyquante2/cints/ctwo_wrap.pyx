cimport cints
from cpython cimport array

def ERI(a,b,c,d):
    if d.contracted:
        return sum(cd*ERI(pd,c,a,b) for (cd,pd) in d)
    return cints.coulomb_repulsion(
    	   a.origin[0],a.origin[1],a.origin[2],a.norm,
	   a.powers[0],a.powers[1],a.powers[2],a.exponent,
    	   b.origin[0],b.origin[1],b.origin[2],b.norm,
	   b.powers[0],b.powers[1],b.powers[2],b.exponent,
    	   c.origin[0],c.origin[1],c.origin[2],c.norm,
	   c.powers[0],c.powers[1],c.powers[2],c.exponent,
    	   d.origin[0],d.origin[1],d.origin[2],d.norm,
	   d.powers[0],d.powers[1],d.powers[2],d.exponent)

