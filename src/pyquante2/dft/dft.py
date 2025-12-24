import numpy as np
from pyquante2.dft.functionals import xs,cvwn5

# Maybe move these to the functionals module and import from there?
xname = dict(lda=xs,xs=xs,svwn=xs)
cname = dict(lda=cvwn5,svwn=cvwn5,xs=None)

def get_xc(grid,D,**kwargs):
    xcname = kwargs.get('xcname','lda')
    # Does not work on either gradient corrected functionals or spin-polarized functionals yet.

    xfunc = xname[xcname]
    cfunc = cname[xcname]
        
    rho = grid.getdens(D)
    fx,dfxa = xfunc(rho)
    if cfunc:
        fc,dfca,dfcb = cfunc(rho,rho)
    else:
        fc=dfca=dfcb=0
        
    w = grid.points[:,3]
    Vxc = np.einsum('g,g,gI,gJ->IJ',w,dfxa+dfca,grid.bfamps,grid.bfamps)
    # The fx comes from either the up or the down spin, whereas the fc comes from
    #  both (which is why x is called with either one, and c is called with both
    Exc = np.dot(w,2*fx+fc)
    return Exc,Vxc
    
