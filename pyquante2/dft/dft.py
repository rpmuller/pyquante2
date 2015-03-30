import numpy as np
from pyquante2.dft.functionals import xs#,cvwn

# Maybe move these to the functionals module and import from there?
xname = dict(lda=xs,xs=xs)
#cname = dict(lda=cvwn)

def get_xc(grid,D,**kwargs):
    xcname = kwargs.get('xcname','lda')
    # Does not work on either gradient corrected functionals or spin-polarized functionals yet.

    xfunctional = xname[xcname]
    #cfunctional = cname[xcname]
        
    rho = grid.getdens(D)
    fx,dfx = xfunctional(rho)
    #fc,dfc = cfunctional(rho)
    # Is it faster to avoid the slice operation on points[:,3] to get the weights, and just
    #  hack the einsum command further?
    w = grid.points[:,3]
    Vx = np.einsum('g,g,gI,gJ->IJ',w,dfx,grid.bfamps,grid.bfamps)
    Ex = np.dot(w,fx)
    return Ex,Vx
    
