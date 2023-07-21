cimport chgp
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
        return chgp.contr_hrr(len(acoefs),a.origin[0],a.origin[1],a.origin[2],anorms.data.as_doubles,
                    a.powers[0],a.powers[1],a.powers[2],aexps.data.as_doubles,acoefs.data.as_doubles,
                    len(bcoefs),b.origin[0],b.origin[1],b.origin[2],bnorms.data.as_doubles,
                    b.powers[0],b.powers[1],b.powers[2],bexps.data.as_doubles,bcoefs.data.as_doubles,
                    len(ccoefs),c.origin[0],c.origin[1],c.origin[2],cnorms.data.as_doubles,
                    c.powers[0],c.powers[1],c.powers[2],cexps.data.as_doubles,ccoefs.data.as_doubles,
                    len(dcoefs),d.origin[0],d.origin[1],d.origin[2],dnorms.data.as_doubles,
                    d.powers[0],d.powers[1],d.powers[2],dexps.data.as_doubles,dcoefs.data.as_doubles)
    if d.contracted:
        return sum(cd*ERI(pd,c,a,b) for (cd,pd) in d)
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
def vrr(double xa,double ya,double za,double norma,int la,int ma,int na,double alphaa,
        double xb,double yb,double zb,double normb,double alphab,
        double xc,double yc,double zc,double normc,int lc,int mc,int nc,double alphac,
        double xd,double yd,double zd,double normd,double alphad,int m):
    return chgp.vrr(xa,ya,za,norma,la,ma,na,alphaa,
                    xb,yb,zb,normb,alphab,
                    xc,yc,zc,normc,lc,mc,nc,alphac,
                    xd,yd,zd,normd,alphad,m)

def vrr_nonrecursive(double xa,double ya,double za,double norma,int la,int ma,int na,double alphaa,
        double xb,double yb,double zb,double normb,double alphab,
        double xc,double yc,double zc,double normc,int lc,int mc,int nc,double alphac,
        double xd,double yd,double zd,double normd,double alphad,int m):
        return chgp.vrr_nonrecursive(xa,ya,za,norma,la,ma,na,alphaa,
				     xb,yb,zb,normb,alphab,
				     xc,yc,zc,normc,lc,mc,nc,alphac,
				     xd,yd,zd,normd,alphad,m)
def vrr_recursive(double xa,double ya,double za,double norma,int la,int ma,int na,double alphaa,
        double xb,double yb,double zb,double normb,double alphab,
        double xc,double yc,double zc,double normc,int lc,int mc,int nc,double alphac,
        double xd,double yd,double zd,double normd,double alphad,int m):
        return chgp.vrr_recursive(xa,ya,za,norma,la,ma,na,alphaa,
				  xb,yb,zb,normb,alphab,
				  xc,yc,zc,normc,lc,mc,nc,alphac,
				  xd,yd,zd,normd,alphad,m)

