/*************************************************************************
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **************************************************************************/

#ifdef _MSC_VER
double lgamma(double);
#endif

static double contr_hrr(int lena, double xa, double ya, double za, double *anorms,
		 int la, int ma, int na, double *aexps, double *acoefs,
		 int lenb, double xb, double yb, double zb, double *bnorms,
		 int lb, int mb, int nb, double *bexps, double *bcoefs,
		 int lenc, double xc, double yc, double zc, double *cnorms,
		 int lc, int mc, int nc, double *cexps, double *ccoefs,
		 int lend, double xd, double yd, double zd, double *dnorms,
		 int ld, int md, int nd, double *dexps, double *dcoefs);

static double contr_vrr(int lena, double xa, double ya, double za,
			double *anorms, int la, int ma, int na,
			double *aexps, double *acoefs,
			int lenb, double xb, double yb, double zb,
			double *bnorms, double *bexps, double *bcoefs,
			int lenc, double xc, double yc, double zc,
			double *cnorms, int lc, int mc, int nc,
			double *cexps, double *ccoefs,
			int lend, double xd, double yd, double zd,
			double *dnorms, double *dexps, double *dcoef);

static double hrr(double xa, double ya, double za, double norma,
	   int la, int ma, int na, double alphaa,
	   double xb, double yb, double zb, double normb,
	   int lb, int mb, int nb, double alphab,
	   double xc, double yc, double zc, double normc,
	   int lc, int mc, int nc, double alphac,
	   double xd, double yd, double zd, double normd,
	   int ld, int md, int nd, double alphad);

static double vrr_recursive(double xa, double ya, double za, double norma,
	   int la, int ma, int na, double alphaa,
	   double xb, double yb, double zb, double normb, double alphab,
	   double xc, double yc, double zc, double normc,
	   int lc, int mc, int nc, double alphac,
	   double xd, double yd, double zd, double normd, double alphad,
	   int m);

static double vrr(double xa, double ya, double za, double norma,
	   int la, int ma, int na, double alphaa,
	   double xb, double yb, double zb, double normb, double alphab,
	   double xc, double yc, double zc, double normc,
	   int lc, int mc, int nc, double alphac,
	   double xd, double yd, double zd, double normd, double alphad,
	   int m);

static int iindex(int la, int ma, int na, int lc, int mc, int nc, int m);

static double dist2(double x1, double y1, double z1, 
		    double x2, double y2, double z2);
static double product_center_1D(double alphaa, double xa, 
				double alphab, double xb);
static double Fgamma(double m, double x);
static double gamm_inc(double a, double x);
static void gser(double *gamser, double a, double x, double *gln);
static void gcf(double *gammcf, double a, double x, double *gln);
