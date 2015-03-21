from pyquante2.dft.functionals import xs

xc_function_map = dict(xs=xs)

def get_xc(grid,D,nel,xcname):
    xc_function = xc_function_map(xcname)
    rho = grid.getdens(D)
    
