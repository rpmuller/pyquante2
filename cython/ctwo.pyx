cimport cints
cimport chgp

STUFF = "Hi" # define init??

# Where do normalization constants come in here?
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

def ERI_hgp(a,b,c,d):
    # This should be faster if I can get it to work, but having trouble passing
    # in to double *anorms, etc. 
    # if a.contracted and b.contracted and c.contracted and d.contracted:
    #     acoefs,anorms,aexps = a.cne_lists()
    #     bcoefs,bnorms,bexps = b.cne_lists()
    #     ccoefs,cnorms,cexps = c.cne_lists()
    #     dcoefs,dnorms,dexps = d.cne_lists()
    #     return chgp.contr_hrr(len(acoefs),a.origin[0],a.origin[1],a.origin[2],anorms2,
    #                           a.powers[0],a.powers[1],a.powers[2],aexps2,acoefs2,
    #                           len(bcoefs),b.origin[0],b.origin[1],b.origin[2],bnorms2,
    #                           b.powers[0],b.powers[1],b.powers[2],bexps2,bcoefs2,
    #                           len(ccoefs),c.origin[0],c.origin[1],c.origin[2],cnorms2,
    #                           c.powers[0],c.powers[1],c.powers[2],cexps2,ccoefs2,
    #                           len(dcoefs),d.origin[0],d.origin[1],d.origin[2],dnorms2,
    #                           d.powers[0],d.powers[1],d.powers[2],dexps2,dcoefs2)
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

