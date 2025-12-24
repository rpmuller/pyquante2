"""
Notes on functionals in pyquante2.

  1. I would like for all operations here to be array-scale operations.
     - This means e.g. np.power rather than pow
  2. I would like as few memory copies as possible.
  3. Need to decide how to handle redundant information, i.e.
     - rhoa and rhob when non-spin-polarized
       This might make tracking total density and the zeta (spin polarization)
       worthwhile; the latter could just be zero (or False or None or whatever) 
     - return values from non-spin-polarized calculations.

"""
import numpy as np
def zero_low_density(rho,cut=1e-10):
    rho[rho<cut]=0
    return rho

def xs(rho,alpha=2/3.):
    "Xalpha X functional. alpha is the X-alpha scaling factor"
    fac=-2.25*alpha*np.power(0.75/np.pi,1./3.)
    rho3 = np.power(rho,1./3.)
    fx = fac*rho*rho3
    dfxdna = (4./3.)*fac*rho3
    return fx,dfxdna

def xb88_array(rho,gam,tol=1e-6):
    # Still doesn't work
    rho = zero_low_density(rho)
    rho13 = np.power(rho,1./3.)
    x = np.zeros(rho.shape,dtype=float)
    g = np.zeros(rho.shape,dtype=float)
    dg = np.zeros(rho.shape,dtype=float)
    x[rho>tol] = np.sqrt(gam)/rho13/rho
    g[rho>tol] = b88_g(x[rho>tol])
    dg[rho>tol] = b88_dg(x[rho>tol])
    dfxdrho = (4./3.)*rho13*(g-x*dg)
    dfxdgam = 0.5*dg/np.sqrt(gam)
    fx = rho*rho13*g
    return fx,dfxdrho,dfxdgam

def xb88(rho,gam,tol=1e-10):
    rho = zero_low_density(rho)
    fxs = []
    dfxdrhos = []
    dfxdgams = []

    for na,gama in zip(rho,gam):
        fx = dfxdrho = dfxdgam = 0
        if na > tol:
            rho13 = np.power(na,1./3.)
            x = np.sqrt(gama)/rho13/na
            g = b88_g(x)
            dg = b88_dg(x)
            dfxdrho = (4./3.)*rho13*(g-x*dg)
            dfxdgam = 0.5*dg/np.sqrt(gama)
            fx = na*rho13*g
        fxs.append(fx)
        dfxdrhos.append(dfxdrho)
        dfxdgams.append(dfxdgam)
    return np.array(fxs),np.array(dfxdrhos),np.array(dfxdgams)

def xpbe(rho,gam,tol=1e-10):
    rho = zero_low_density(rho)
    fxs = []
    dfxdrhos = []
    dfxdgams = []

    for na,gama in zip(rho,gam):
        fx = dfxdrho = dfxdgam = 0
        if na > tol:
            kap = 0.804
            mu = 0.449276922095889E-2
            fx0,vx0 = xs(na)
            rho13 = na**(1.0/3.0)
            rho43 = rho13*na
            den = 1+mu*gama/rho43/rho43
            F = 1+kap-kap/den
            fx = fx0*F
            dFdr = -(8./3.)*kap*mu*gama/den/den*na**(-11./3.)
            dfxdrho = vx0*F+fx0*dFdr
            dFdg = -kap*mu/rho43/rho43/den/den
            dfxdgam = fx0*dFdg
        fxs.append(fx)
        dfxdrhos.append(dfxdrho)
        dfxdgams.append(dfxdgam)
    return np.array(fxs),np.array(dfxdrhos),np.array(dfxdgams)

def cvwn5(rhoa,rhob,tol=1e-10):
    rhoa = zero_low_density(rhoa)
    rhob = zero_low_density(rhob)

    ecs = []
    vcrhoas = []
    vcrhobs = []
    for na,nb in zip(rhoa,rhob):
        rho = na+nb
        ec = vcrhoa = vcrhob = 0
        if rho>tol:
            zeta=(na-nb)/rho
            x = pow(3./4./np.pi/rho,1/6.)
            epsp = vwn_epsp(x)
            epsf = vwn_epsf(x)
            g = vwn_g(zeta)
            eps = epsp + g*(epsf-epsp)
            ec = eps*rho

            depsp = vwn_depsp(x)
            depsf = vwn_depsf(x)
            dg = vwn_dg(zeta)
            deps_dx = depsp + g*(depsf-depsp)
            deps_dg = (epsf-epsp)*dg
            vcrhoa = eps - (x/6.)*deps_dx + deps_dg*(1-zeta)
            vcrhob = eps - (x/6.)*deps_dx - deps_dg*(1+zeta)
        ecs.append(ec)
        vcrhoas.append(vcrhoa)
        vcrhobs.append(vcrhob)
    return np.array(ecs),np.array(vcrhoas),np.array(vcrhobs)

def clyp(rhoas,rhobs,gaas,gabs,gbbs,tol=1e-10):
    fcs = []
    fcnas = []
    fcnbs = []
    fcgaas = []
    fcgabs = []
    fcgbbs = []
    for na,nb,gaa,gab,gbb in zip(rhoas,rhobs,gaas,gabs,gbbs):
        fc,fcna,fcnb,fcgaa,fcgab,fcgbb = clyp_point(na,nb,gaa,gab,gbb,tol)
        fcs.append(fc)
        fcnas.append(fcnbs)
        fcnbs.append(fcnb)
        fcgaas.append(fcgaa)
        fcgabs.append(fcgab)
        fcgbbs.append(fcgbb)
    return np.array(fcs),np.array(fcnas),np.array(fcnbs),np.array(fcgaas),np.array(fcgabs),np.array(fcgbbs)

def clyp_point(rhoa,rhob,gamaa,gamab,gambb,tol=1e-10):
    # Modified and corrected by AEM in June 2006.
    a = 0.04918  # Parameters from the LYP papers
    b = 0.132
    c = 0.2533
    d = 0.349
    rho = rhoa+rhob
    fc=fcrhoa=fcrhob=fcgamaa=fcgamab=fcgambb=0
    assert rhoa >= 0.0
    assert rhob >= 0.0
    if rho > tol:
        rhom3 = np.power(rho,-1./3.)
        w = np.exp(-c*rhom3)/(1+d*rhom3)*np.power(rho,-11./3.)
        dl = c*rhom3+d*rhom3/(1+d*rhom3)

        fcgamaa = -a*b*w*((1./9.)*rhoa*rhob*(1-3*dl-(dl-11)*rhoa/rho)-rhob*rhob)
        fcgamab = -a*b*w*((1./9.)*rhoa*rhob*(47-7*dl)-(4./3.)*rho*rho)
        fcgambb = -a*b*w*((1./9.)*rhoa*rhob*(1-3*dl-(dl-11)*rhob/rho)-rhoa*rhoa)

        fc = -4*a/(1+d*rhom3)*rhoa*rhob/rho \
             -np.power(2,11./3.)*0.3*np.power(3*np.pi*np.pi,2./3.)*a*b*w \
             *rhoa*rhob*(np.power(rhoa,8./3.)+np.power(rhob,8./3.)) \
             + fcgamaa*gamaa + fcgamab*gamab + fcgambb*gambb

        dw = -(1./3.)*np.power(rho,-4./3.)*w*(11*np.power(rho,1./3.)-c-d/(1+d*rhom3))
        ddl = (1./3.)*(d*d*np.power(rho,-5./3.)/np.power(1+d*rhom3,2)-dl/rho)

        d2f_dradgaa = dw/w*fcgamaa - a*b*w*(
            (1./9.)*rhob*(1-3*dl-(dl-11)*rhoa/rho)
            -(1./9.)*rhoa*rhob*((3+rhoa/rho)*ddl+(dl-11)*rhob/rho/rho))
        d2f_dradgbb = dw/w*fcgambb - a*b*w*(
            (1./9.)*rhob*(1-3*dl-(dl-11)*rhob/rho)
            -(1./9.)*rhoa*rhob*((3+rhob/rho)*ddl-(dl-11)*rhob/rho/rho)
            -2*rhoa)
        d2f_dradgab = dw/w*fcgamab-a*b*w*(
            (1./9)*rhob*(47-7*dl)-(7./9.)*rhoa*rhob*ddl-(8./3.)*rho)

        d2f_drbdgaa = dw/w*fcgamaa - a*b*w*(
            (1./9.)*rhoa*(1-3*dl-(dl-11)*rhoa/rho)
            -(1./9.)*rhoa*rhob*((3+rhoa/rho)*ddl-(dl-11)*rhoa/rho/rho)
            -2*rhob)
        d2f_drbdgbb = dw/w*fcgambb - a*b*w*(
            (1./9.)*rhoa*(1-3*dl-(dl-11)*rhob/rho)
            -(1./9.)*rhoa*rhob*((3+rhob/rho)*ddl+(dl-11)*rhoa/rho/rho))
        d2f_drbdgab = dw/w*fcgamab-a*b*w*(
            (1./9)*rhoa*(47-7*dl)-(7./9.)*rhoa*rhob*ddl-(8./3.)*rho)

        fcrhoa = fcrhob = 0
        if rhoa > tol:
            fcrhoa = -4*a/(1+d*rhom3)*rhoa*rhob/rho*(
                (1./3.)*d*np.power(rho,-4./3.)/(1+d*rhom3)+1/rhoa-1/rho)\
                -np.power(2,11./3.)*0.3*np.power(3*np.pi*np.pi,2./3.)*a*b*(
                dw*rhoa*rhob*(np.power(rhoa,8./3.)+np.power(rhob,8./3.))
                +w*rhob*((11./3.)*np.power(rhoa,8./3.)+np.power(rhob,8./3.))) \
                +d2f_dradgaa*gamaa + d2f_dradgbb*gambb + d2f_dradgab*gamab

        if rhob > tol:
            fcrhob = -4*a/(1+d*rhom3)*rhoa*rhob/rho*(
                (1./3.)*d*np.power(rho,-4./3.)/(1+d*rhom3)+1/rhob-1/rho)\
                -np.power(2,11./3.)*0.3*np.power(3*np.pi*np.pi,2./3.)*a*b*(
                dw*rhoa*rhob*(np.power(rhob,8./3.)+np.power(rhoa,8./3.))
                +w*rhoa*((11./3.)*np.power(rhob,8./3.)+np.power(rhoa,8./3.))) \
                +d2f_drbdgaa*gamaa + d2f_drbdgbb*gambb + d2f_drbdgab*gamab
    return fc,fcrhoa,fcrhob,fcgamaa,fcgamab,fcgambb

def cpbe(na,nb,ga,gab,gb):
    "PBE Correlation Functional"
    npts = len(na)
    ec = np.zeros(npts,'d')
    vca = np.zeros(npts,'d')
    vcb = np.zeros(npts,'d')
    vcga = np.zeros(npts,'d')
    vcgab = np.zeros(npts,'d')
    vcgb = np.zeros(npts,'d')
    for i in range(npts):
        ec[i],vca[i],vcb[i],vcga[i],vcgab[i],vcgb[i] = \
                            cpbe_point(na[i],nb[i],ga[i],gab[i],gb[i])
    return ec,vca,vcb,vcga,vcgab,vcgb

def cpbe_point(rhoa,rhob,gama,gamb,gamab,tol=1e-10):
    rho = rhoa+rhob
    ec = vca = vcb = vcgama = vcgamb = vcgamab = 0
    gam = 0.031091
    ohm = 0.046644
    bet = 0.066725
    if rho > tol:
        Rs = np.power(3./(4.*np.pi*rho),1./3.)
        Zeta = (rhoa-rhob)/rho
        Kf = np.power(3*np.pi*np.pi*rho,1./3.)
        Ks = np.sqrt(4*Kf/np.pi)
        Phi = 0.5*(np.power(1+Zeta,2./3.) + np.power(1-Zeta,2./3.))
        Phi3 = Phi*Phi*Phi
        gradrho = np.sqrt(gama+gamb+2.*gamab)
        T = gradrho/(2*Phi*Ks*rho)
        T2 = T*T
        T4 = T2*T2

        eps,vc0a,vc0b = cpbe_lsd(rhoa,rhob)

        expo = (np.exp(-eps/(gam*Phi3))-1.)
        A = bet/gam/expo
        N = T2+A*T4
        D = 1.+A*T2+A*A*T4
        H = gam*Phi3*np.log(1.+(bet/gam)*N/D)
        ec = rho*(eps+H)

        # Derivative stuff
        dZ_drhoa = (1.-Zeta)/rho
        dZ_drhob = -(1.+Zeta)/rho

        dPhi_dZ = np.power(1.+Zeta,-1./3.)/3.-np.power(1.-Zeta,-1./3.)/3.
        dPhi_drhoa = dPhi_dZ*dZ_drhoa
        dPhi_drhob = dPhi_dZ*dZ_drhob
        
        dKs_drho = Ks/(6*rho)
        
        dT_dPhi = -T/Phi
        dT_dKs = -T/Ks
        dT_drhoa = -T/rho + dT_dPhi*dPhi_drhoa + dT_dKs*dKs_drho
        dT_drhob = -T/rho + dT_dPhi*dPhi_drhob + dT_dKs*dKs_drho

        dA_dPhi = -A/expo*np.exp(-eps/(gam*Phi3))*(3*eps/(gam*Phi3*Phi))
        dA_deps = -A/expo*np.exp(-eps/(gam*Phi3))*(-1/(gam*Phi3))
        deps_drhoa = (vc0a-eps)/rho
        deps_drhob = (vc0b-eps)/rho
        dA_drhoa = dA_dPhi*dPhi_drhoa + dA_deps*deps_drhoa
        dA_drhob = dA_dPhi*dPhi_drhob + dA_deps*deps_drhoa

        dN_dT = 2*T+4*A*T2*T
        dD_dT = 2*A*T + 4*A*A*T*T2
        dN_dA = T4
        dD_dA = T2+2*A*T4

        dH_dPhi = 3*H/Phi
        dH_dT = bet*Phi3/(1.+bet/gam*N/D)*(D*dN_dT-N*dD_dT)/D/D
            
        dH_dA = bet*Phi3/(1.+bet/gam*N/D)*(D*dN_dA-N*dD_dA)/D/D
        
        dH_drhoa = dH_dPhi*dPhi_drhoa + dH_dT*dT_drhoa + dH_dA*dA_drhoa
        dH_drhob = dH_dPhi*dPhi_drhob + dH_dT*dT_drhob + dH_dA*dA_drhob
        
        vca = vc0a + H + rho*dH_drhoa
        vcb = vc0b + H + rho*dH_drhob
    # Havent done the dE_dgamma derives yet
    return ec,vca,vcb,vcgama,vcgamab,vcgamb

def vwn_xx(x,b,c): return x*x+b*x+c
def vwn_epsp(x): return vwn_eps(x,0.0310907,-0.10498,3.72744,12.9352)
#def vwn_epsf(x): return vwn_eps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_epsf(x): return vwn_eps(x,0.01554535,-0.32500,7.06042,18.0578)

def vwn_eps(x,a,x0,b,c):
    Q = np.sqrt(4*c-b*b)
    eps = a*(np.log(x*x/vwn_xx(x,b,c))
             - b*(x0/vwn_xx(x0,b,c))*np.log(np.power(x-x0,2)/vwn_xx(x,b,c))
             + (2*b/Q)*(1-(x0*(2*x0+b)/vwn_xx(x0,b,c))) * np.arctan(Q/(2*x+b)))
    #eps = a*(np.log(x*x/vwn_xx(x,b,c)) + (2*b/Q)*np.arctan(Q/(2*x+b))
    #         - (b*x0/vwn_xx(x0,b,c))*np.log(np.power(x-x0,2)/vwn_xx(x,b,c))
    #         + (2*(b+2*x0)/Q)*np.arctan(Q/(2*x+b)))
    return eps

def vwn_eps0(x,a,x0,b,c):
    def X(x): return x*x+b*x+c
    Q = np.sqrt(4*c-b*b)
    eps = a*(np.log(x*x/X(x)) + (2*b/Q)*np.arctan(Q/(2*x+b))
             - (b*x0/X(x0))*np.log(np.power(x-x0,2)/X(x))
             + (2*(b+2*x0)/Q)*np.arctan(Q/(2*x+b)))
    return eps

def vwn_depsp(x): return vwn_deps(x,0.0310907,-0.10498,3.72744,12.9352)
#def vwn_depsf(x): return vwn_deps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_depsf(x): return vwn_deps(x,0.01554535,-0.32500,7.06042,18.0578)
def vwn_deps(x,a,x0,b,c):
    q = np.sqrt(4*c-b*b)
    deps = a*(2/x - (2*x+b)/vwn_xx(x,b,c)
              - 4*b/(np.power(2*x+b,2)+q*q) - (b*x0/vwn_xx(x0,b,c))
              * (2/(x-x0)-(2*x+b)/vwn_xx(x,b,c)-4*(2*x0+b)/(np.power(2*x+b,2)+q*q)))
    return deps

def vwn_g(z): return 1.125*(np.power(1+z,4./3.)+np.power(1-z,4./3.)-2)
def vwn_dg(z): return 1.5*(np.power(1+z,1./3.)-np.power(1-z,1./3.))

def b88_g(x,b=0.0042):
    return -1.5*np.power(3./4./np.pi,1./3.)-b*x*x/(1.+6.*b*x*np.arcsinh(x))

def b88_dg(x,b=0.0042):
    num = 6*b*b*x*x*(x/np.sqrt(x*x+1)-np.arcsinh(x))-2*b*x
    denom = np.power(1+6*b*x*np.arcsinh(x),2)
    return num/denom

def cpbe_lsd(rhoa,rhob):
    # Not quite VWN. AEM: It's usually called PW correlation
    # LSD terms
    # Note that this routine gives out ec, not fc.
    # If you rather have fc, use pw instead
    rho = rhoa+rhob
    Rs = np.power(3./(4.*np.pi*rho),1./3.)
    Zeta = (rhoa-rhob)/rho
    thrd = 1./3.     # thrd*=various multiples of 1/3
    thrd4 = 4*thrd
    ggam=0.5198420997897463295344212145565 # gam= 2^(4/3)-2
    fzz=8./(9.*ggam) # fzz=f''(0)= 8/(9*gam)
    rtrs = np.sqrt(Rs)
    eu,eurs = pbe_gcor(0.0310907,0.21370,7.5957,
                       3.5876,1.6382,0.49294,rtrs)
    ep,eprs = pbe_gcor(0.01554535,0.20548,14.1189,
                       6.1977,3.3662,0.62517,rtrs)
    alfm,alfrsm = pbe_gcor(0.0168869,0.11125,10.357,
                           3.6231,0.88026,0.49671,rtrs)
    alfc = -alfm
    z4 = Zeta**4
    f=(np.power(1.+Zeta,thrd4)+np.power(1.-Zeta,thrd4)-2.)/ggam
    eps = eu*(1.-f*z4)+ep*f*z4-alfm*f*(1.-z4)/fzz

    ecrs = eurs*(1.-f*z4)+eprs*f*z4-alfrsm*f*(1.-z4)/fzz
    fz = thrd4*(np.power(1.+Zeta,thrd)-np.power(1.-Zeta,thrd))/ggam
    eczet = 4.*(Zeta**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu-(1.-z4)*alfm/fzz)
    comm = eps -Rs*ecrs/3.-Zeta*eczet
    vca = comm + eczet
    vcb = comm - eczet
    return eps,vca,vcb
    
def pbe_gcor(a,a1,b1,b2,b3,b4,rtrs):
#      subroutine gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
# slimmed down version of gcor used in pw91 routines, to interpolate
# lsd correlation energy, as given by (10) of
# j. p. perdew and y. wang, phys. rev. b {\bf 45}, 13244 (1992).
# k. burke, may 11, 1996.
#      implicit real*8 (a-h,o-z)
      q0 = -2.*a*(1.+a1*rtrs*rtrs)
      q1 = 2.*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = np.log(1.+1./q1)
      gg = q0*q2
      q3 = a*(b1/rtrs+2.*b2+rtrs*(3.*b3+4.*b4*rtrs))
      ggrs = -2.*a*a1*q2-q0*q3/(q1*(1.+q1))
      return gg,ggrs

if __name__ == '__main__':
    import pylab as pyl
    x = np.linspace(1e-12,1)
    pyl.plot(x,vwn_eps(x,0.0310907,-0.10498,3.72744,12.9352))
    pyl.plot(x,vwn_eps0(x,0.0310907,-0.10498,3.72744,12.9352))
    pyl.show()
