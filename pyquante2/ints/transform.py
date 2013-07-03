import numpy as np
def transform(ints,c):
    return np.einsum('aI,bJ,cK,dL,abcd->IJKL',c,c,c,c,ints)

def transform_mp2(ints,c,nocc):
    return np.einsum('aI,bJ,cK,dL,abcd->IJKL',c[:,:nocc],c,c[:,:nocc],c,ints)
    
