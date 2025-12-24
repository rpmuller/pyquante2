cdef extern from "crys.h":

    double contr_coulomb(int lena,double *aexps,double *acoefs,double *anorms,
		     double xa,double ya,double za,int la,int ma,int na,
		     int lenb,double *bexps,double *bcoefs,double *bnorms,
		     double xb,double yb,double zb,int lb,int mb,int nb,
		     int lenc,double *cexps,double *ccoefs,double *cnorms,
		     double xc,double yc,double zc,int lc,int mc,int nc,
		     int lend,double *dexps,double *dcoefs,double *dnorms,
		     double xd,double yd,double zd,int ld,int md,int nd)

    double coulomb_repulsion(double xa,double ya,double za,double norma,
			 int la,int ma,int na,double alphaa,
			 double xb,double yb,double zb,double normb,
			 int lb,int mb,int nb,double alphab,
			 double xc,double yc,double zc,double normc,
			 int lc,int mc,int nc,double alphac,
			 double xd,double yd,double zd,double normd,
			 int ld,int md,int nd,double alphad)

