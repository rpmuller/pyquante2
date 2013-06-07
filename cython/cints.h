/*************************************************************************
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **************************************************************************/
double fB(int i, int l1, int l2, double px, double ax, double bx, 
	  int r, double g);
double Bfunc(int i, int r, double g);
double contr_coulomb(int ia, double *aexps, double *acoefs, double *anorms,
			    double xa, double ya, double za, int la, int ma, int na, 
			    int ib, double *bexps, double *bcoefs, double *bnorms,
			    double xb, double yb, double zb, int lb, int mb, int nb, 
			    int ic, double *cexps, double *ccoefs, double *cnorms,
			    double xc, double yc, double zc, int lc, int mc, int nc, 
			    int id, double *dexps, double *dcoefs, double *dnorms,
			    double xd, double yd, double zd, int ld, int md, int nd);

double coulomb_repulsion(double xa, double ya, double za, double norma,
				int la, int ma, int na, double alphaa,
				double xb, double yb, double zb, double normb,
				int lb, int mb, int nb, double alphab,
				double xc, double yc, double zc, double normc,
				int lc, int mc, int nc, double alphac,
				double xd, double yd, double zd, double normd,
				int ld, int md, int nd, double alphad);

double *B_array(int l1, int l2, int l3, int l4, double p, double a,
		double b, double q, double c, double d,
		double g1, double g2, double delta);

double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
		     int l3, int l4, double Px, double Ax, double Bx,
		     double Qx, double Cx, double Dx, double gamma1,
		     double gamma2, double delta);
double kinetic(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb);
double overlap(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb);
double overlap_1D(int l1, int l2, double PAx,
		  double PBx, double gamma);
double nuclear_attraction(double x1, double y1, double z1, double norm1,
				 int l1, int m1, int n1, double alpha1,
				 double x2, double y2, double z2, double norm2,
				 int l2, int m2, int n2, double alpha2,
				 double x3, double y3, double z3);
double A_term(int i, int r, int u, int l1, int l2,
		     double PAx, double PBx, double CPx, double gamma);
double *A_array(int l1, int l2, double PA, double PB,
		       double CP, double g);

int fact(int n);
int fact2(int n);
double dist2(double x1, double y1, double z1, 
		    double x2, double y2, double z2);
double dist(double x1, double y1, double z1, 
		   double x2, double y2, double z2);
double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb);
int binomial(int a, int b);

double Fgamma(double m, double x);
double gamm_inc(double a, double x);

int ijkl2intindex(int i, int j, int k, int l);

int fact_ratio2(int a, int b);

double product_center_1D(double alphaa, double xa, 
			 double alphab, double xb);

double three_center_1D(double xi, int ai, double alphai,
			      double xj, int aj, double alphaj,
			      double xk, int ak, double alphak);

/* Routines from Numerical Recipes */
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
