cimport cints
cimport chgp
import ctypes
from cpython cimport array

STUFF = "Hi" # define init??

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

# Having a contracted-only version of the function doesn't help.
def cERI_hgp(a,b,c,d):
    cdef array.array acoefs2,anorms2,aexps2
    cdef array.array bcoefs2,bnorms2,bexps2
    cdef array.array ccoefs2,cnorms2,cexps2
    cdef array.array dcoefs2,dnorms2,dexps2
    acoefs,anorms,aexps = a.cne_list()	    
    bcoefs,bnorms,bexps = b.cne_list()	    
    ccoefs,cnorms,cexps = c.cne_list()	    
    dcoefs,dnorms,dexps = d.cne_list()	    
    acoefs2 = array.array('d',acoefs)
    bcoefs2 = array.array('d',bcoefs)
    ccoefs2 = array.array('d',ccoefs)
    dcoefs2 = array.array('d',dcoefs)
    anorms2 = array.array('d',anorms)
    bnorms2 = array.array('d',bnorms)
    cnorms2 = array.array('d',cnorms)
    dnorms2 = array.array('d',dnorms)
    aexps2 = array.array('d',aexps)
    bexps2 = array.array('d',bexps)
    cexps2 = array.array('d',cexps)
    dexps2 = array.array('d',dexps)
    return chgp.contr_hrr(len(acoefs),a.origin[0],a.origin[1],a.origin[2],anorms2.data.as_doubles,
                a.powers[0],a.powers[1],a.powers[2],aexps2.data.as_doubles,acoefs2.data.as_doubles,
                len(bcoefs),b.origin[0],b.origin[1],b.origin[2],bnorms2.data.as_doubles,
                b.powers[0],b.powers[1],b.powers[2],bexps2.data.as_doubles,bcoefs2.data.as_doubles,
                len(ccoefs),c.origin[0],c.origin[1],c.origin[2],cnorms2.data.as_doubles,
                c.powers[0],c.powers[1],c.powers[2],cexps2.data.as_doubles,ccoefs2.data.as_doubles,
                len(dcoefs),d.origin[0],d.origin[1],d.origin[2],dnorms2.data.as_doubles,
                d.powers[0],d.powers[1],d.powers[2],dexps2.data.as_doubles,dcoefs2.data.as_doubles)

def ERI_hgp(a,b,c,d):
    # This should be faster if I can get it to work, but having trouble passing
    # in to double *anorms, etc.
    cdef array.array acoefs2,anorms2,aexps2
    cdef array.array bcoefs2,bnorms2,bexps2
    cdef array.array ccoefs2,cnorms2,cexps2
    cdef array.array dcoefs2,dnorms2,dexps2
    if a.contracted and b.contracted and c.contracted and d.contracted:
        acoefs,anorms,aexps = a.cne_list()	    
        bcoefs,bnorms,bexps = b.cne_list()	    
        ccoefs,cnorms,cexps = c.cne_list()	    
        dcoefs,dnorms,dexps = d.cne_list()	    
        acoefs2 = array.array('d',acoefs)
        bcoefs2 = array.array('d',bcoefs)
        ccoefs2 = array.array('d',ccoefs)
        dcoefs2 = array.array('d',dcoefs)
        anorms2 = array.array('d',anorms)
        bnorms2 = array.array('d',bnorms)
        cnorms2 = array.array('d',cnorms)
        dnorms2 = array.array('d',dnorms)
        aexps2 = array.array('d',aexps)
        bexps2 = array.array('d',bexps)
        cexps2 = array.array('d',cexps)
        dexps2 = array.array('d',dexps)
        return chgp.contr_hrr(len(acoefs),a.origin[0],a.origin[1],a.origin[2],anorms2.data.as_doubles,
                    a.powers[0],a.powers[1],a.powers[2],aexps2.data.as_doubles,acoefs2.data.as_doubles,
                    len(bcoefs),b.origin[0],b.origin[1],b.origin[2],bnorms2.data.as_doubles,
                    b.powers[0],b.powers[1],b.powers[2],bexps2.data.as_doubles,bcoefs2.data.as_doubles,
                    len(ccoefs),c.origin[0],c.origin[1],c.origin[2],cnorms2.data.as_doubles,
                    c.powers[0],c.powers[1],c.powers[2],cexps2.data.as_doubles,ccoefs2.data.as_doubles,
                    len(dcoefs),d.origin[0],d.origin[1],d.origin[2],dnorms2.data.as_doubles,
                    d.powers[0],d.powers[1],d.powers[2],dexps2.data.as_doubles,dcoefs2.data.as_doubles)
    if d.contracted:
        return sum(cd*ERI_hgp(pd,c,a,b) for (cd,pd) in d)
    return chgp.hrr(
        a.origin[0],a.origin[1],a.origin[2],a.norm,
        a.powers[0],a.powers[1],a.powers[2],a.exponent,
	b.origin[0],b.origin[1],b.origin[2],b.norm,
	b.powers[0],b.powers[1],b.powers[2],b.exponent,
	c.origin[0],c.origin[1],c.origin[2],c.norm,
	c.powers[0],c.powers[1],c.powers[2],c.exponent,
	d.origin[0],d.origin[1],d.origin[2],d.norm,
	d.powers[0],d.powers[1],d.powers[2],d.exponent)

# The following are only for debugging and can be deleted after ERI_hgp works:
def vrr(xa,ya,za,norma,la,ma,na,alphaa,
        xb,yb,zb,normb,alphab,
        xc,yc,zc,normc,lc,mc,nc,alphac,
        xd,yd,zd,normd,alphad,m):
    return chgp.vrr(xa,ya,za,norma,la,ma,na,alphaa,
                    xb,yb,zb,normb,alphab,
                    xc,yc,zc,normc,lc,mc,nc,alphac,
                    xd,yd,zd,normd,alphad,m)

def vrr_nonrecursive(
	xa,ya,za,norma,la,ma,na,alphaa,
	xb,yb,zb,normb,alphab,
	xc,yc,zc,normc,lc,mc,nc,alphac,
	xd,yd,zd,normd,alphad,m):
        return chgp.vrr_nonrecursive(xa,ya,za,norma,la,ma,na,alphaa,
				     xb,yb,zb,normb,alphab,
				     xc,yc,zc,normc,lc,mc,nc,alphac,
				     xd,yd,zd,normd,alphad,m)
def vrr_recursive(
	xa,ya,za,norma,la,ma,na,alphaa,
	xb,yb,zb,normb,alphab,
	xc,yc,zc,normc,lc,mc,nc,alphac,
	xd,yd,zd,normd,alphad,m):
        return chgp.vrr_recursive(xa,ya,za,norma,la,ma,na,alphaa,
				  xb,yb,zb,normb,alphab,
				  xc,yc,zc,normc,lc,mc,nc,alphac,
				  xd,yd,zd,normd,alphad,m)

