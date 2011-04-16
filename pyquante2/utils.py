"""
utils.py - Simple utilility funtions used in pyquante2.

"""

from math import sqrt,factorial,log,exp

def fact2(n):
    """
    fact2(n) - n!!, double factorial of n
    >>> fact2(0)
    1
    >>> fact2(1)
    1
    >>> fact2(3)
    3
    >>> fact2(8)
    384
    >>> fact2(-1)
    1
    """
    return reduce(int.__mul__,xrange(n,0,-2),1)

def dist2(a,b):
    """
    Square of the cartesian distance.
    
    >>> dist2((0,0,0),(1,0,0))
    1
    >>> dist2((0,0,0),(0,2,0))
    4
    >>> dist2((4,0,0),(0,0,0))
    16
    >>> dist2((0,5,0),(0,0,0))
    25
    >>> dist2((2,2,2),(0,0,0))
    12
    """
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2

def dist(a,b):
    """
    The Cartesian distance
    >>> dist((0,0,0),(0,0,0))
    0.0
    >>> dist((0,0,0),(1,0,0))
    1.0
    """
    return sqrt(dist2(a,b))

def binomial(n,k):
    """
    Binomial coefficient
    >>> binomial(5,2)
    10
    >>> binomial(10,5)
    252
    """
    if n==k: return 1
    assert n>k, "Attempting to call binomial(%d,%d)" % (n,k)
    return factorial(n)/(factorial(k)*factorial(n-k))

def Fgamma(m,x):
    """
    Incomplete gamma function
    >>> Fgamma(0,0)
    1.0000000000000016
    """
    SMALL=1e-15
    x = max(abs(x),SMALL)
    val = gamm_inc(m+0.5,x)
    return 0.5*pow(x,-m-0.5)*val;

def gammln(x):
    """
    Numerical recipes, section 6.1
    >>> gammln(1)
    0.0
    """
    cof = [76.18009172947146,-86.50532032941677,
           24.01409824083091,-1.231739572450155,
           0.1208650973866179e-2,-0.5395239384953e-5]
    y=x
    tmp=x+5.5
    tmp = tmp - (x+0.5)*log(tmp)
    ser=1.000000000190015 # don't you just love these numbers?!
    for j in xrange(6):
        y = y+1
        ser = ser+cof[j]/y
    return -tmp+log(2.5066282746310005*ser/x);

def gamm_inc(a,x):
    """
    Incomple gamma function \gamma; computed from NumRec routine gammp.
    """
    gammap,gln = gammp(a,x)
    return exp(gln)*gammap
    
def gammp(a,x):
    "Returns the incomplete gamma function P(a;x). NumRec sect 6.2."
    assert (x > 0 and a >= 0), "Invalid arguments in routine gammp"

    if x < (a+1.0): #Use the series representation
        gamser,gln = _gser(a,x)
        return gamser,gln
    #Use the continued fraction representation
    gammcf,gln = _gcf(a,x)
    return 1.0-gammcf ,gln  #and take its complement.

def _gser(a,x):
    "Series representation of Gamma. NumRec sect 6.1."
    ITMAX=100
    EPS=3.e-7

    gln=gammln(a)
    assert(x>=0),'x < 0 in gser'
    if x == 0 : return 0,gln

    ap = a
    delt = sum = 1./a
    for i in xrange(ITMAX):
        ap=ap+1.
        delt=delt*x/ap
        sum=sum+delt
        if abs(delt) < abs(sum)*EPS: break
    else:
        print 'a too large, ITMAX too small in gser'
    gamser=sum*exp(-x+a*log(x)-gln)
    return gamser,gln

def _gcf(a,x):
    "Continued fraction representation of Gamma. NumRec sect 6.1"
    ITMAX=100
    EPS=3.e-7
    FPMIN=1.e-30

    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    for i in xrange(1,ITMAX+1):
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if abs(d) < FPMIN: d=FPMIN
        c=b+an/c
        if abs(c) < FPMIN: c=FPMIN
        d=1./d
        delt=d*c
        h=h*delt
        if abs(delt-1.) < EPS: break
    else:
        print 'a too large, ITMAX too small in gcf'
    gammcf=exp(-x+a*log(x)-gln)*h
    return gammcf,gln
    
    
def erf(z): 
    if z==0.0: return 0.0
    else: return gammp(0.5,z**2)[0]

if __name__ == '__main__':
    import doctest
    doctest.testmod()

