/**********************************************************************
 * cints.c  C implementation of simple math functions in pyutil.
 *
 * The equations herein are based upon
 * 'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
 * S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
 * [THO paper].
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **********************************************************************/

#include "cints.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define SMALL 0.00000001

double fB(int i, int l1, int l2, double px, double ax, double bx, 
		 int r, double g){
  return binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,r,g);
}

double Bfunc(int i, int r, double g){
  return fact_ratio2(i,r)*pow(4*g,r-i);
}

double contr_coulomb(int lena, double *aexps, double *acoefs,
			    double *anorms, double xa, double ya, double za,
			    int la, int ma, int na, 
			    int lenb, double *bexps, double *bcoefs,
			    double *bnorms, double xb, double yb, double zb,
			    int lb, int mb, int nb, 
			    int lenc, double *cexps, double *ccoefs,
			    double *cnorms, double xc, double yc, double zc,
			    int lc, int mc, int nc, 
			    int lend, double *dexps, double *dcoefs,
			    double *dnorms, double xd, double yd, double zd,
			    int ld, int md, int nd){

  int i,j,k,l;
  double Jij = 0.,incr=0.;

  for (i=0; i<lena; i++)
    for (j=0; j<lenb; j++)
      for (k=0; k<lenc; k++)
	for (l=0; l<lend; l++){
	  incr = coulomb_repulsion(xa,ya,za,anorms[i],la,ma,na,aexps[i],
			      xb,yb,zb,bnorms[j],lb,mb,nb,bexps[j],
			      xc,yc,zc,cnorms[k],lc,mc,nc,cexps[k],
			      xd,yd,zd,dnorms[l],ld,md,nd,dexps[l]);
	  
	  Jij += acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*incr;
	}
  return Jij;
}

double coulomb_repulsion(double xa, double ya, double za, double norma,
				int la, int ma, int na, double alphaa,
				double xb, double yb, double zb, double normb,
				int lb, int mb, int nb, double alphab,
				double xc, double yc, double zc, double normc,
				int lc, int mc, int nc, double alphac,
				double xd, double yd, double zd, double normd,
				int ld, int md, int nd, double alphad){

  double rab2, rcd2,rpq2,xp,yp,zp,xq,yq,zq,gamma1,gamma2,delta,sum;
  double *Bx, *By, *Bz;
  int I,J,K;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  rcd2 = dist2(xc,yc,zc,xd,yd,zd);
  xp = product_center_1D(alphaa,xa,alphab,xb);
  yp = product_center_1D(alphaa,ya,alphab,yb);
  zp = product_center_1D(alphaa,za,alphab,zb);
  xq = product_center_1D(alphac,xc,alphad,xd);
  yq = product_center_1D(alphac,yc,alphad,yd);
  zq = product_center_1D(alphac,zc,alphad,zd);
  rpq2 = dist2(xp,yp,zp,xq,yq,zq);
  gamma1 = alphaa+alphab;
  gamma2 = alphac+alphad;
  delta = (1./gamma1+1./gamma2)/4.;

  Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta);
  By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta);
  Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta);

  sum = 0.;
  for (I=0; I<la+lb+lc+ld+1;I++)
    for (J=0; J<ma+mb+mc+md+1;J++)
      for (K=0; K<na+nb+nc+nd+1;K++)
	sum += Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta);

  free(Bx);
  free(By);
  free(Bz);  
  
  return 2.*pow(M_PI,2.5)/(gamma1*gamma2*sqrt(gamma1+gamma2))
    *exp(-alphaa*alphab*rab2/gamma1) 
    *exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd;
}

double *B_array(int l1, int l2, int l3, int l4, double p, double a,
		double b, double q, double c, double d,
		double g1, double g2, double delta){
  int Imax,i1,i2,r1,r2,u,I,i;
  double *B;
  Imax = l1+l2+l3+l4+1;
  B = (double *)malloc(Imax*sizeof(double));
  for (i=0; i<Imax; i++) B[i] = 0.;

  for (i1=0; i1<l1+l2+1; i1++)
    for (i2=0; i2<l3+l4+1; i2++)
      for (r1=0; r1<i1/2+1; r1++)
	for (r2=0; r2<i2/2+1; r2++)
	  for (u=0; u<(i1+i2)/2-r1-r2+1; u++){
	    I = i1+i2-2*(r1+r2)-u;
	    B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
				 p,a,b,q,c,d,g1,g2,delta);
	  }

  return B;
}

double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
	      int l3, int l4, double Px, double Ax, double Bx,
	      double Qx, double Cx, double Dx, double gamma1,
	      double gamma2, double delta){
  /* THO eq. 2.22 */
  return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)
    *pow(-1,i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)
    *pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u)
    *pow(Qx-Px,i1+i2-2*(r1+r2)-2*u)
    /pow(delta,i1+i2-2*(r1+r2)-u);
}


double kinetic(double alpha1, int l1, int m1, int n1,
	       double xa, double ya, double za,
	       double alpha2, int l2, int m2, int n2,
	       double xb, double yb, double zb){

  double term0,term1,term2;
  term0 = alpha2*(2*(l2+m2+n2)+3)*
    overlap(alpha1,l1,m1,n1,xa,ya,za,
		   alpha2,l2,m2,n2,xb,yb,zb);
  term1 = -2*pow(alpha2,2)*
    (overlap(alpha1,l1,m1,n1,xa,ya,za,
		    alpha2,l2+2,m2,n2,xb,yb,zb)
     + overlap(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2+2,n2,xb,yb,zb)
     + overlap(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2,n2+2,xb,yb,zb));
  term2 = -0.5*(l2*(l2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2-2,m2,n2,xb,yb,zb) +
		m2*(m2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2-2,n2,xb,yb,zb) +
		n2*(n2-1)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2,n2-2,xb,yb,zb));
  return term0+term1+term2;
}

double overlap(double alpha1, int l1, int m1, int n1,
		      double xa, double ya, double za,
		      double alpha2, int l2, int m2, int n2,
		      double xb, double yb, double zb){
  /*Taken from THO eq. 2.12*/
  double rab2,gamma,xp,yp,zp,pre,wx,wy,wz;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  gamma = alpha1+alpha2;
  xp = product_center_1D(alpha1,xa,alpha2,xb);
  yp = product_center_1D(alpha1,ya,alpha2,yb);
  zp = product_center_1D(alpha1,za,alpha2,zb);

  pre = pow(M_PI/gamma,1.5)*exp(-alpha1*alpha2*rab2/gamma);

  wx = overlap_1D(l1,l2,xp-xa,xp-xb,gamma);
  wy = overlap_1D(m1,m2,yp-ya,yp-yb,gamma);
  wz = overlap_1D(n1,n2,zp-za,zp-zb,gamma);
  return pre*wx*wy*wz;
}

double overlap_1D(int l1, int l2, double PAx,
		  double PBx, double gamma){
  /*Taken from THO eq. 2.12*/
  int i;
  double sum;
  sum = 0.;
  for (i=0; i<(1+floor(0.5*(l1+l2))); i++)
    sum += binomial_prefactor(2*i,l1,l2,PAx,PBx)* 
      fact2(2*i-1)/pow(2*gamma,i);
  return sum;
}
    
double nuclear_attraction(double x1, double y1, double z1, double norm1,
				 int l1, int m1, int n1, double alpha1,
				 double x2, double y2, double z2, double norm2,
				 int l2, int m2, int n2, double alpha2,
				 double x3, double y3, double z3){
  int I,J,K;
  double gamma,xp,yp,zp,sum,rab2,rcp2;
  double *Ax,*Ay,*Az;

  gamma = alpha1+alpha2;

  xp = product_center_1D(alpha1,x1,alpha2,x2);
  yp = product_center_1D(alpha1,y1,alpha2,y2);
  zp = product_center_1D(alpha1,z1,alpha2,z2);

  rab2 = dist2(x1,y1,z1,x2,y2,z2);
  rcp2 = dist2(x3,y3,z3,xp,yp,zp);

  Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma);
  Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma);
  Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma);

  sum = 0.;
  for (I=0; I<l1+l2+1; I++)
    for (J=0; J<m1+m2+1; J++)
      for (K=0; K<n1+n2+1; K++)
	sum += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma);

  free(Ax);
  free(Ay);
  free(Az);
  return -norm1*norm2*
    2*M_PI/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum;
}
    
double A_term(int i, int r, int u, int l1, int l2,
		     double PAx, double PBx, double CPx, double gamma){
  /* THO eq. 2.18 */
  return pow(-1,i)*binomial_prefactor(i,l1,l2,PAx,PBx)*
    pow(-1,u)*fact(i)*pow(CPx,i-2*r-2*u)*
    pow(0.25/gamma,r+u)/fact(r)/fact(u)/fact(i-2*r-2*u);
}

double *A_array(int l1, int l2, double PA, double PB,
		double CP, double g){
  /* THO eq. 2.18 and 3.1 */
  int Imax,i,r,u,I;
  double *A;

  Imax = l1+l2+1;
  A = (double *)malloc(Imax*sizeof(double));
  for (i=0; i<Imax; i++) A[i] = 0.;
  for (i=0; i<Imax; i++)
    for (r=0; r<floor(i/2)+1;r++)
      for (u=0; u<floor((i-2*r)/2.)+1; u++){
	I = i-2*r-u;
	A[I] += A_term(i,r,u,l1,l2,PA,PB,CP,g);
      }
  return A;
}


int fact(int n){
  if (n <= 1) return 1;
  return n*fact(n-1);
}

int fact2(int n){ /* double factorial function = 1*3*5*...*n */
  if (n <= 1) return 1;
  return n*fact2(n-2);
}

double dist2(double x1, double y1, double z1,
		    double x2, double y2, double z2){
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}
double dist(double x1, double y1, double z1,
		   double x2, double y2, double z2){
  return sqrt(dist2(x1,y1,z1,x2,y2,z2));
}

double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb){
  int t;
  double sum=0.;
  for (t=0; t<s+1; t++)
    if ((s-ia <= t) && (t <= ib)) 
      sum += binomial(ia,s-t)*binomial(ib,t)*pow(xpa,ia-s+t)*pow(xpb,ib-t);
  return sum;
} 

int binomial(int a, int b){return fact(a)/(fact(b)*fact(a-b));}

double Fgamma(double m, double x){
  double val;
  if (fabs(x) < SMALL) x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5*pow(x,-m-0.5)*val; 
}

double gamm_inc(double a, double x){ /* Taken from NR routine gammap */
  double gamser,gammcf,gln;
  
  assert (x >= 0.);
  assert (a > 0.);
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return exp(gln)*gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return exp(gln)*(1.0-gammcf);
  }
}
 
void gser(double *gamser, double a, double x, double *gln){
  int n;
  double sum,del,ap;

  *gln=lgamma(a);
  if (x <= 0.0) {
    assert(x>=0.);
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser");
    return;
  }
}
 
void gcf(double *gammcf, double a, double x, double *gln){
  int i;
  double an,b,c,d,del,h;
  
  *gln=lgamma(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  assert(i<=ITMAX);
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

int ijkl2intindex(int i, int j, int k, int l){
  int tmp,ij,kl;
  if (i<j) return ijkl2intindex(j,i,k,l);
  if (k<l) return ijkl2intindex(i,j,l,k);
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij<kl){
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

int fact_ratio2(int a, int b){ return fact(a)/fact(b)/fact(a-2*b); }

double product_center_1D(double alphaa, double xa, 
			 double alphab, double xb){
  return (alphaa*xa+alphab*xb)/(alphaa+alphab);
}

double three_center_1D(double xi, int ai, double alphai,
			      double xj, int aj, double alphaj,
			      double xk, int ak, double alphak){

  double gamma, dx, px, xpi,xpj,xpk,intgl;
  int q,r,s,n;
  
  gamma = alphai+alphaj+alphak;
  dx = exp(-alphai*alphaj*pow(xi-xj,2)/gamma) *
    exp(-alphai*alphak*pow(xi-xk,2)/gamma) *
    exp(-alphaj*alphak*pow(xj-xk,2)/gamma);
  px = (alphai*xi+alphaj*xj+alphak*xk)/gamma;
    
  xpi = px-xi;
  xpj = px-xj;
  xpk = px-xk;
  intgl = 0;
  for (q=0; q<ai+1; q++){
    for (r=0; r<aj+1; r++){
      for (s=0; s<ak+1; s++){
	n = (q+r+s)/2;
	if ((q+r+s)%2 == 0) {
	  intgl += binomial(ai,q)*binomial(aj,r)*binomial(ak,s)*
	    pow(xpi,ai-q)*pow(xpj,aj-r)*pow(xpk,ak-s)*
	    fact2(2*n-1)/pow(2*gamma,n)*sqrt(M_PI/gamma);
	}
      }
    }
  }
  return dx*intgl;
}


/* work is the work space for the various exponents, contraction */
/*  coefficients, etc., used by the contracted code. Decided to  */
/*  allocate this all at once rather than doing mallocs/frees. */
/*  Only speeds things up a little, but every little bit counts, I guess. */
#define MAX_PRIMS_PER_CONT (10)
double work[12*MAX_PRIMS_PER_CONT];

#undef ITMAX
#undef EPS
#undef FPMIN
