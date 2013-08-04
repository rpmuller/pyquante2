import numpy as np
from itertools import product
from functools import reduce

def mp2(ints,orbs,orbe,nocc,nvirt,verbose=False):
    moints = ints.transform_mp2(orbs,nocc)
    Emp2 = 0
    for a,b in product(range(nocc),repeat=2):
        Eab = 0
        for r,s in product(range(nocc,nocc+nvirt),repeat=2):
            arbs,asbr = moints[a,r,b,s],moints[a,s,b,r]
            Eab += arbs*(2*arbs-asbr)/ \
                   (orbe[a]+orbe[b]-orbe[r]-orbe[s])
        if verbose: print("MP2 pair energy for %d,%d: %f" % (a,b,Eab))
        Emp2 += Eab
    if verbose: print("MP2 energy: %f" % Eab)
    return Emp2

if __name__ == '__main__': test_mp2()
