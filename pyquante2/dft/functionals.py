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

def xb88(rho,gam):
    rho13 = pow(rho,1./3.)
    x = np.sqrt(gam)/rho13/rho 
    g = b88_g(x)
    dg = b88_dg(x)
    dfxdrho = (4./3.)*rho13*(g-x*dg)
    dfxdgam = 0.5*dg/np.sqrt(gam)
    fx = rho*rho13*g
    return fx,dfxdrho,dfxdgam

def xpbe(rho,gam):
    kap = 0.804
    mu = 0.449276922095889E-2
    ex0,vx0 = xs(rho)
    rho13 = rho**(1.0/3.0)
    rho43 = rho13*rho
    den = 1.E0+mu*gam/rho43/rho43
    F = 1+kap-kap/den
    ex = ex0*F
    dFdr = -(8./3.)*kap*mu*gam/den/den*rho**(-11./3.)
    vxrho = vx0*F+ex0*dFdr
    dFdg = -kap*mu/rho43/rho43/den/den
    vxgam = ex0*dFdg
    return ex,vxrho,vxgam

def cvwn_vector(rhoa,rhob):
    rho = rhoa+rhob
    zeta=(rhoa-rhob)/rho
    x = np.power(3./4./np.pi/rho,1./6.)
    epsp = vwn_epsp(x)
    depsp = vwn_depsp(x)
    g = vwn_g(zeta)
    epsf = vwn_epsf(x)
    eps = epsp + g*(epsf-epsp)
    ec = eps*rho
    depsf = vwn_depsf(x)
    dg = vwn_dg(zeta)
    deps_dx = depsp + g*(depsf-depsp)
    deps_dg = (epsf-epsp)*dg
    vcrhoa = eps - (x/6.)*deps_dx + deps_dg*(1-zeta)
    vcrhob = eps - (x/6.)*deps_dx - deps_dg*(1+zeta)
    return ec,vcrhoa,vcrhob

def cvwn_pointwise(rhoa,rhob,tol=1e-10):
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
            x = pow(3./4./np.pi/rho,1./6.)
            epsp = vwn_epsp(x)
            depsp = vwn_depsp(x)
            g = vwn_g(zeta)
            epsf = vwn_epsf(x)
            eps = epsp + g*(epsf-epsp)
            ec = eps*rho
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

def cvwn(rhoa,rhob): return cvwn_pointwise(rhoa,rhob)


def vwn_xx(x,b,c): return x*x+b*x+c
def vwn_epsp(x): return vwn_eps(x,0.0310907,-0.10498,3.72744,12.9352)
def vwn_epsf(x): return vwn_eps(x,0.01554535,-0.32500,7.06042,13.0045)
def vwn_eps(x,a,x0,b,c):
    q = np.sqrt(4*c-b*b)
    eps = a*(np.log(x*x/vwn_xx(x,b,c))
             - b*(x0/vwn_xx(x0,b,c))*np.log(np.power(x-x0,2)/vwn_xx(x,b,c))
             + (2*b/q)*(1-(x0*(2*x0+b)/vwn_xx(x0,b,c))) * np.arctan(q/(2*x+b)))
    return eps

def vwn_depsp(x): return vwn_deps(x,0.0310907,-0.10498,3.72744,12.9352)
def vwn_depsf(x): return vwn_deps(x,0.01554535,-0.32500,7.06042,13.0045)
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

if __name__ == '__main__':
    import pylab as pyl
    z = np.linspace(0,1)
    pyl.plot(z,vwn_g(z))
    pyl.plot(z,vwn_dg(z))
    pyl.show()
