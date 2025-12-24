cimport crys
from cpython cimport array

def ERI(a,b,c,d):
    cdef array.array acoefs,anorms,aexps
    cdef array.array bcoefs,bnorms,bexps
    cdef array.array ccoefs,cnorms,cexps
    cdef array.array dcoefs,dnorms,dexps
    if a.contracted and b.contracted and c.contracted and d.contracted:
        acoefs,anorms,aexps = a.cne_list()
        bcoefs,bnorms,bexps = b.cne_list()
        ccoefs,cnorms,cexps = c.cne_list()
        dcoefs,dnorms,dexps = d.cne_list()
        return crys.contr_coulomb(
	    len(acoefs),aexps.data.as_doubles,acoefs.data.as_doubles,anorms.data.as_doubles,
            a.origin[0],a.origin[1],a.origin[2],a.powers[0],a.powers[1],a.powers[2],
	    len(bcoefs),bexps.data.as_doubles,bcoefs.data.as_doubles,bnorms.data.as_doubles,
            b.origin[0],b.origin[1],b.origin[2],b.powers[0],b.powers[1],b.powers[2],
	    len(ccoefs),cexps.data.as_doubles,ccoefs.data.as_doubles,cnorms.data.as_doubles,
            c.origin[0],c.origin[1],c.origin[2],c.powers[0],c.powers[1],c.powers[2],
	    len(dcoefs),dexps.data.as_doubles,dcoefs.data.as_doubles,dnorms.data.as_doubles,
            d.origin[0],d.origin[1],d.origin[2],d.powers[0],d.powers[1],d.powers[2])
    if d.contracted:
        return sum(cd*ERI(pd,c,a,b) for (cd,pd) in d)
    return crys.coulomb_repulsion(
        a.origin[0],a.origin[1],a.origin[2],a.norm,a.powers[0],a.powers[1],a.powers[2],a.exponent,
        b.origin[0],b.origin[1],b.origin[2],b.norm,b.powers[0],b.powers[1],b.powers[2],b.exponent,
        c.origin[0],c.origin[1],c.origin[2],c.norm,c.powers[0],c.powers[1],c.powers[2],c.exponent,
        d.origin[0],d.origin[1],d.origin[2],d.norm,d.powers[0],d.powers[1],d.powers[2],d.exponent)

