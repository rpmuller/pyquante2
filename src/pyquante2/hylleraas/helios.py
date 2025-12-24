#!/usr/bin/env python
"""\
Helios.py: Solve two-electron atoms using Pekeris techniques
Copyright 1999, 2000 Richard P. Muller, David R. Kent IV, and
   William A. Goddard, III

Special thanks to Edward Montgomery, Michael Barnett, and
Robert Forrey, without whom this work would have been impossible.


Test values (from Pekeris) for verification.
             Singlet states:
n  Matrix Size     H-                He
5    22        0.52763068142   2.90368898612
9    95        0.52775001651   2.90372338908
10   125       0.52775061025   2.90372387862
11   161       0.52775085979   2.90372411115
12   203       0.52775093560   2.90372422832
13   252       0.52775097384   2.90372429041
16   444       0.52775100630   2.90372435622
19   715       0.52775101339   2.90372437081
22   1078      0.52775101536   2.90372437476
             Triplet states:
n  Matrix Size    He
11   125       2.17522097961
14   252       2.17522925889
17   444       2.17522937679

>>> np.isclose(two_electron_solve(1,5,0),-0.52763068141531577)
True
>>> np.isclose(two_electron_solve(2,5,0),-2.9036889861207293)
True
>>> np.isclose(two_electron_solve(1,9,0),-0.52775001651538511)
True
>>> np.isclose(two_electron_solve(2,9,0),-2.9037233890716729)
True
>>> np.isclose(two_electron_solve(2,11,1),-2.1752209796141289)
True
"""
import numpy as np
from pyquante2.utils import geigh

def kval(l,m,n,spin):
    # This is from Pekeris, except that the last term of the singlet
    # expression corrects a type in Pekeris' paper
    w = l+m+n
    lm = l+m
    if spin == 0:
        k = w*(w+2)*(2*w+5)/24. + (1-pow(-1,w))/16. + lm*lm/4. +\
            (1-pow(-1,lm))/8. + l+ 1 + lm/2.
    elif spin == 1:
        # The following was caught by Ed Montgomery, 1/16/00. Note
        # the different sign from the previous expression.
        #k = w*(w+2)*(2*w-1)/24. + (1-pow(-1,w))/16. + l*(m+n) + m
        k = w*(w+2)*(2*w-1)/24. - (1-pow(-1,w))/16. + l*(m+n) + m
    else:
        print("kval: Error -- unknown spin")
        sys.exit()
    # k should now be an integer
    k = int(k)
    return k

def korder(wmax,spin):
    # Return a list of tuples with (l,m,n) where:
    # l    Exponent of s term
    # m    Exponent of t term
    # n    Exponent of u term
    klist = []
    if spin == 0:
        for w in range(wmax):
            for l in range(w+1):
                for m in range(w+1):
                    n = w-l-m
                    if n>=0 and l<=m:
                        klist.append((l,m,n))
    elif spin == 1:
        for w in range(wmax):
            for n in range(w):
                for m in range(w+1):
                    l = w-n-m
                    if 0<=l<m:
                        klist.append((l,m,n))
    else:
        print("korder: ERROR improper spin state")
        sys.exit()
    return klist

def hterm(l,m,n,l2,m2,n2):
    # Obtain the Pekeris Hamiltonian:
    # Hc = ScE
    # H = a*Z + b
    # S = c
    delta = l2-l, m2-m, n2-n
    if delta == (2,0,0):
        x = 4*(l+1)*(l+2)
        a = -x
        b = 0
        c = x*(1+m+n)
    elif delta == (0,2,0):
        x = 4*(m+1)*(m+2)
        a = -x
        b = 0
        c = x*(1+l+n)
    elif delta == (1,1,0):
        x = 4*(l+1)*(m+1)
        a = -2*x
        b = x
        c = x*(2+l+m)
    elif delta == (1,0,1):
        x = 2*(l+1)*(n+1)
        a = -2*x
        b = x
        c = x*(2+2*m+n)
    elif delta == (0,1,1):
        x = 2*(m+1)*(n+1)
        a = -2*x
        b = x
        c = x*(2+2*l+n)
    elif delta == (0,0,2):
        x = (n+1)*(n+2)
        a = 0
        b = x
        c = 0
    elif delta == (1,0,0):
        x = l+1
        a = 4*x*(4*l+4*m+2*n+7)
        b = x*(-8*m-4*n-6)
        c = -2*x*((m+n)*(4*m+12*l)+n*n+12*l+18*m+15*n+14)
    elif delta == (0,1,0):
        x = m+1
        a = 4*x*(4*l+4*m+2*n+7)
        b = x*(-8*l-4*n-6)
        c = -2*x*((l+n)*(4*l+12*m)+n*n+12*m+18*l+15*n+14)
    elif delta == (0,0,1):
        x = 4*(n+1)
        a = x*(2*l+2*m+2)
        b = x*(-l-m-n-2)
        c = -x*(-l*l-m*m+4*l*m+2*l*n+2*n*m+3*l+3*m+2*n+2)
    elif delta == (0,2,-1):
        x = 4*(m+1)*(m+2)*n
        a = 0
        b = 0
        c = x
    elif delta == (2,0,-1):
        x = 4*(l+1)*(l+2)*n
        a = 0
        b = 0
        c = x
    elif delta == (-1,0,2):
        x = 2*l*(n+1)*(n+2)
        a = 0
        b = 0
        c = x
    elif delta == (0,-1,2):
        x = 2*m*(n+1)*(n+2)
        a = 0
        b = 0
        c = x
    elif delta == (0,0,0):
        a = -4*((l+m)*(6*l+6*m+4*n+12)-4*l*m+4*n+8)
        b = 4*(2*l+1)*(2*m+1)+4*(2*n+1)*(l+m+1)+6*n*n+6*n+2
        c = 4*((l+m)*(10*l*m+10*m*n+10*l*n+10*l+10*m+18*n+4*n*n+16) +\
               l*m*(4-12*n)+8+12*n+4*n*n)
    elif delta == (-1,1,0):
        x = 4*l*(m+1)
        a = -2*x
        b = x
        c = x*(1+l+m)
    elif delta == (1,-1,0):
        x = 4*(l+1)*m
        a = -2*x
        b = x
        c = x*(1+l+m)
    elif delta == (-1,0,1):
        x = 2*l*(n+1)
        a = -2*x
        b = x
        c = x*(2*m-4*l-n)
    elif delta == (0,-1,1):
        x = 2*m*(n+1)
        a = -2*x
        b = x
        c = x*(2*l-4*m-n)
    elif delta == (1,0,-1):
        x = 2*(l+1)*n
        a = -2*x
        b = x
        c = x*(2*m-4*l-n-3)
    elif delta == (0,1,-1):
        x = 2*(m+1)*n
        a = -2*x
        b = x 
        c = x*(2*l-4*m-n-3)
    elif delta == (-1,0,0):
        x = 2*l
        a = x*(8*l+8*m+4*n+6)
        b = -x*(4*m+2*n+3)
        c = -x*((m+n+1)*(12*l+4*m+2)+n+n*n)
    elif delta == (0,-1,0):
        x = 2*m
        a = x*(8*l+8*m+4*n+6)
        b = -x*(4*l+2*n+3)
        c = -x*((l+n+1)*(12*m+4*l+2)+n+n*n)
    elif delta == (0,0,-1):
        x = 4*n
        a = x*(2*l+2*m+2)
        b = -x*(l+m+n+1)
        c = -x*((l+m)*(1+2*n-l-m)+6*l*m+2*n)
    elif delta == (1,0,-2):
        x = 2*n*(n-1)*(l+1)
        a = 0
        b = 0
        c = x
    elif delta == (0,1,-2):
        x = 2*n*(n-1)*(m+1)
        a = 0
        b = 0
        c = x
    elif delta == (-2,0,1):
        x = 4*l*(l-1)*(n+1)
        a = 0
        b = 0
        c = x
    elif delta == (0,-2,1):
        x = 4*m*(m-1)*(n+1)
        a = 0
        b = 0
        c = x
    elif delta == (-2,0,0):
        x = 4*l*(l-1)
        a = -x
        b = 0
        c = x*(1+m+n)
    elif delta == (0,-2,0):
        x = 4*m*(m-1)
        a = -x
        b = 0
        c = x*(1+l+n)
    elif delta == (0,0,-2):
        x = n*(n-1)
        a = 0
        b = x
        c = 0
    elif delta == (-1,-1,0):
        x = 4*l*m
        a = -2*x
        b = x
        c = x*(l+m)
    elif delta == (-1,0,-1):
        x = 2*l*n
        a = -2*x
        b = x
        c = x*(2*m+n+1)
    elif delta == (0,-1,-1):
        x = 2*m*n
        a = -2*x
        b = x
        c = x*(2*l+n+1)
    else:
        a = 0.
        b = 0.
        c = 0.
    return (a,b,c)

def pekeris(Z,wmax,spin):
    # Return Pekeris H and S of order n with nuclear charge Z
    # write H = a*Z + b

    klist = korder(wmax,spin)
    N = len(klist)
    H = np.zeros((N,N),dtype=float)
    S = np.zeros((N,N),dtype=float)
    for index1 in range(N):
        l,m,n = klist[index1]
        k = kval(l,m,n,spin)
        for index2 in range(N):
            l2,m2,n2 = klist[index2]
            k2 = kval(l2,m2,n2,spin)
            i = k-1
            j = k2-1
            a,b,c = hterm(l,m,n,l2,m2,n2)
            if l == m and l2 == m2:
                a = 0.5*a
                b = 0.5*b
                c = 0.5*c
            elif l == m or l2 == m2:
                pass #do nothing here
            elif spin == 1:
                a2,b2,c2 = hterm(m,l,n,l2,m2,n2)
                a = a - a2
                b = b - b2
                c = c - c2
            elif spin == 0:
                a2,b2,c2 = hterm(m,l,n,l2,m2,n2)
                a = a + a2
                b = b + b2
                c = c + c2
            else:
                print("pekeris: ERROR should not be here")
                sys.exit()
            H[i,j] = a*Z + b
            S[i,j] = c
                
    return (H,S)

def transform(A,B):
    # Similarity transformation: returns (B+)AB
    C = matrixmultiply(A,B)
    return matrixmultiply(transpose(B),C)

def inv_sqrt(M):
    # Returns the inverse square root of a matrix
    E,V = eigenvectors(M)
    n = len(E)
    M = zeros((n,n),Float)
    for i in range(n):
        M[i,i] = 1./sqrt(E[i])
    return transform(M,V)

def two_electron_solve(atomic_number,maximum_order,spin):
    Z = atomic_number
    H,S = pekeris(Z,maximum_order,spin)
    E,V = geigh(H,S)
    epsilon = E[0]
    E2 = [-Ei**2 for Ei in E]
    #print("Energy (h) for order %d: %15.12f %15.12f" % (len(E),E2[0],epsilon))
    return E2[0]

# Main program starts here:
if __name__ == '__main__':
    import doctest; doctest.testmod()
