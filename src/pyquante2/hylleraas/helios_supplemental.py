#!/usr/bin/env python

# This file contains supplemental (and deprecated) routines from the
# helios program on two-electron atoms.


def hterm_1962(l,m,n,l2,m2,n2):
    # Obtain the Pekeris Hamiltonian:
    # Hc = ScE
    # H = a*Z + b
    # S = c
    # This routine is a variant of the normal hterm one, from the
    # 1962 paper (PR 126,1470 (1962))
    delta = l2-l,m2-m,n2-n
    a = 0.
    b = 0.
    c = 0.
    if delta == (0,-1,0):
        c = 4*(l*(l+1)-(m-1)*(2*n+1)+n*(n-1))
    elif delta == (-1,0,0):
        c = 4*(m*(m+1)-(l-1)*(2*n+1)+n*(n-1))
    elif delta == (1,-2,0):
        c = 4*(l+1)*(l+1)
    elif delta == (-2,1,0):
        c = 4*(m+1)*(m+1)
    elif delta == (1,-1,-1):
        c = 4*(l+1)*(l+1)
    elif delta == (-1,1,-1):
        c = 4*(m+1)*(m+1)
    elif delta == (0,0,-1):
        c = 4*(l*l + m*m - (n-1)*(l+m+1))
    elif delta == (1,0,-2):
        c = 2*(l+1)*(l+1)
    elif delta == (0,1,-2):
        c = 2*(m+1)*(m+1)
    elif delta == (-2,0,1):
        c = 8*(n+1)*(n+1)
    elif delta == (0,-2,1):
        c = 8*(n+1)*(n+1)
    elif delta == (-1,-1,0):
        a = 8
        b = -4
        c = -4*(l+m)
    elif delta == (0,-2,0):
        a = 4
        c = -4*(l+n+1)
    elif delta == (-2,0,0):
        a = 4
        b = -4*(m+n+1)
    elif delta == (0,-1,-1):
        a = 4
        b = -2
        c = -2*(2*l+n+1)
    elif delta == (-1,0,-1):
        a = 4
        b = -2
        c = -2*(2*m+n+1)
    elif delta == (-2,-1,0):
        c = 2
    elif delta == (-1,-2,0):
        c = 2
    elif delta == (-1,-1,-1):
        c = 2
    elif delta == (0,-2,-1):
        c = 1
    elif delta == (-2,0,-1):
        c = 1
    elif delta == (-1,0,-2):
        c = 0.5
    elif delta == (0,-1,-2):
        c = 0.5
    elif delta == (0,0,-2):
        b = -1
    return (a,b,c)


def hterm_cox_sutcliffe(l,m,n,l2,m2,n2):
    # Obtain the Pekeris Hamiltonian:
    # Hc = ScE
    # H = a*Z + b
    # S = c
    K = 1
    M = 2
    delta = l2-l, m2-m, n2-n
    a = 0.
    b = 0.
    c = 0.
    x = 0.
    if delta == (2,1,0):
        c = 2.*(l+2)*(l+1)*(K-1)*(m+1)
    elif delta == (2,0,1):
        c = (l+2)*(l+1)(K*M-M-2)*(n+1)/M
    elif delta == (1,2,0):
        c = 2*(m+2)*(m+1)*(K-1)*(l+1)
    elif delta == (1,0,2):
        c = (n+2)*(n+1)*(K*M-M-2)*(l+1)/(2*M)
    elif delta == (0,2,1):
        c = (m+2)*(m+1)*(K*M-M-2)*(n+1)/M
    elif delta == (0,1,2):
        c = (n+2)*(n+1)*(KM-M-2)*(m+1)/(2*M)
    elif delta == (1,1,1):
        c = 2*(K*M-M-2)*(l+1)*(n+1)*(m+1)/M
    elif delta == (2,0,0):
        x = -(l+2)*(l+1)/M
        a = x*(-4*M)
        c = x*(2*n*M + M + 4*n + 2 + 3*K*M + 4*K*m*M + 2*K*n*M)
    elif delta == (0,2,0):
        x = -(m+2)*(m+1)/M
        a = x*(-4*M)
        c = x*(2*K*n*M + 2*n*M + 3*K*M + 4*n + 2 + M + 4*K*l*M)
    elif delta == (0,0,2):
        x = -(n+2)*(n+1)/M
        b = x
        c = x*(-2*l-M-2*m-2+K*M-m*M-l*M+K*m*M+K*l*M)
    elif delta == (1,1,0):
        x = 2*(l+1)*(m+1)/M
        a = x*(-4*M)
        b = x*2*M
        c = x*(4*K*m*M-2*m*M-2*l*M+4*K*l*M+4*n+2
               +2*K*n*M-5*M-2*n*M+9*K*M)
    elif delta == (1,0,1):
        x = 2*(l+1)*(n+1)
        a = x*(-2*M)
        b = x*M
        c = x*(K*n*M-2*l*M+2*K*l*M-4*l+4*K*M-4*m+2*K*m*M-6-2*M)
    elif delta == (0,1,1):
        x = 2*(n+1)*(m+1)
        a = x*(-2*M)
        b = x*M
        c = x*(2*K*m*M-4*m-2*m*M+K*n*M-6-4*l+2*K*l*M-2*M+4*K*M)
    elif delta == (1,0,0):
        x = (l+1)/M
        a = x*(-28*M-16*m*M-8*n*M-16*l*M)
        b = x*(6*M + 8*m*M + 4*n*M)
        c = x*(14+12*m+12*l+16*l*n+8*m*l+9*M+2*n*n+4*m*M+26*n
               +15*n*M+19*K*M+3*K*n*n*M+16*m*n-4*m*m*M+12*K*m*m*M
               +16*K*l*m*M+8*K*l*n*M+8*K*m*n*M+32*m*M+15*K*n*M
               +12*K*l*M-4*m*m+16*l*n*M+8*m*l*M+12*l*M-n*n*M)
    elif delta == (0,1,0):
        x = (m+1)/M
        a = x*(-28*M-16*m*M-8*n*M-16*l*M)
        b = x*(6*M+4*n*M+8*l*M)
        c = x*(14+12*m+12*l+16*l*n+8*m*l+9*M+2*n*n+12*m*M+26*n
               -4*l*l+15*n*M+19*K*M+3*K*n*n*M+16*m*n+16*K*l*m*M
               +16*m*n*M+8*K*l*n*M+8*m*n*M+12*K*m*M+15*K*n*M
               +32*K*l*M+12*K*l*l*M+8*m*l*M+8*l*M+4*l*M-n*n*M-4*l*l*M)
    elif delta == (0,0,1):
        x = 2*(n+1)/M
        a = x*(-4*M-4*m*M-4*l*M)
        b = x*(4*M+2*m*M+2*n*M+2*l*M)
        c = x*(-6-10*m-10*l-8*m*l+2*m*M-M-m*M-6*l*l+2*n*M+5*K*M
               -5*m*m*M+3*K*m*m*M+4*K*l*m*M+2*m*n*M+2*K*l*n*M
               +2*K*m*n*M+7*K*m*M+2*K*n*M+7*K*l*M-6*m*m
               +3*K*l*l*M+2*l*n*M+4*m*l*M-l*M-4*l*l*M)
    elif delta == (2,-1,0):
        c = 2*m*(l+2)*(l+1)*(K*M-2-M)/M
    elif delta == (-1,2,0):
        c = 2*l*(m+2)*(m+1)*(K*M-2-M)/M
    elif delta == (0,2,-1):
        c = n*(m+2)*(m+1)*(K*M-5*M-2)/M
    elif delta == (2,0,-1):
        c = n*(l+2)*(l+1)*(K*M-5*M-2)/M
    elif delta == (-1,0,2):
        c = l*(n+2)*(n+1)*(K*M-5*M-2)/(2*M)
    elif delta == (0,-1,2):
        c = m*(n+2)*(n+1)*(K*M-5*M-2)/(2*M)
    elif delta == (1,1,-1):
        c = 2*n*(K*M-2-M)*(l+1)*(m+1)/M
    elif delta == (1,-1,1):
        c = 2*m*(K*M-2-M)*(l+1)*(n+1)/M
    elif delta == (-1,1,1):
        c = 2*l*(K*M-2-M)*(n+1)*(m+1)/M
    elif delta == (0,0,0):
        x = 1./M
        a = x*(-32*M-48*l*M-48*m*M-24*m*m*M-24*l*l*M-16*n*M
             -32*l*m*M-16*l*n*M-16*m*n*M)
        b = x*(-16*l*m*M+12*m*M+14*n*M+12*l*M+6*n*n*M+10*M
               +8*m*n*M+8*l*n*M)
        c = x*(4*m*n*n+44*l*n+32*m+32*l+4*n*n+20*l*l+28*n+20*m*m
               +16+16*K*M+32*l*m*n+16*M+44*m*n+16*l*l*m+16*l*m*m
               +48*l*m+24*l*l*n+4*l*n*n+24*m*m*n+32*m*M+32*l*M
               +10*n*n*M+22*l*l*M+30*n*M+22*m*m*M+10*m*n*n*M
               +32*l*M+32*K*m*M+28*m*m*n*M+18*K*l*l*M+6*K*n*n*M
               +18*K*n*M)
    elif delta == (-1,1,0):
        x = -2*l*(m+1)/M
        a = x*(-4*M)
        b = x*(2*M)
        c = x*(4*K*m*M-2*m*M-4*m+4*K*l*M+5*K*M-2-2*n*M
               -2*l*M+2*K*n*M+4*n-3*M-4*l)
    elif delta == (1,-1,0):
        x = -2*m*(l+1)/M
        a = x*(-4*M)
        b = x*(2*M)
        c = x*(4*K*l*M-2*l*M-4*l-2-4*m+2*K*n*M+4*n-3*M
               +5*K*M-2*m*M-2*n*M+4*K*m*M)
    elif delta == (-1,0,1):
        x = -2*l*(n+1)/M
        a = x*(-2*M)
        b = x
        c = x*(K*n*M-2*n*M+2*K*l*M+2*K*M-4*m-2-6*l*M+2*K*m*M-4*l-2*M)
    elif delta == (0,-1,1):
        x = -2*m*(n+1)/M
        a = x*(-2*M)
        b = x
        c = x*(K*n*M-2*n*M+2*K*m*M-4*l+2*K*M-2-6*m*M+2*K*l*M-2*M-4*m)
    elif delta == (1,0,-1):
        x = -2*n*(l+1)/M
        a = x*(-2*M)
        b = x
        c = x*(2*K*l*M-6*l*M-4*l+K*n*M-6-6*M+2*K*m*M+3*K*M-4*m-2*n*M)
    elif delta == (0,1,-1):
        x = -2*n*(m+1)/M
        a = x*(-2*M)
        b = x
        c = x*(2*K*m*M-6*m*M-4*m-6-6*M-4*l+K*n*M+3*K*M+2*K*l*M-2*n*M)
    elif delta == (-1,0,0):
        x = l/M
        a = x*(-12*M-16*m*M-8*n*M-16*l*M)
        b = x*(6*M+4*n*M+8*m*M)
        c = x*(2+4*m+8*m*l+10*n-3*M+12*l-4*m*m*M-4*m*M
               +16*m*n+3*K*n*n*M-n*n*M-4*m*m+12*K*l*M
               +8*K*l*n*M+7*K*M+7*K*n*M+16*l*n*M-n*M
               +12*l*M+8*K*m*n*M+16*l*n+8*m*l*M+16*K*m*M
               +12*K*m*m*M+16*K*l*m*M+2*n*n)
    elif delta == (0,-1,0):
        x = m/M
        a = x*(-12*M-16*m*M-8*n*M-16*l*M)
        b = x*(6*M+4*n*M+8*l*M)
        c = x*(2+12*m+8*m*l+10*n-3*M+4*l-4*l*l+12*m*M
               +16*m*n+3*K*n*n*M-n*n*M+16*K*l*M+8*K*l*n*M
               +7*K*M+7*K*n*M-n*M-4*l*M+8*K*m*n*M+16*l*n+8*m*l*M
               +12*K*m*M+16*K*l*m*M-4*l*l*M+2*n*n+12*K*l*l*M+16*m*n*M)
    elif delta == (0,0,-1):
        x = 2*n/M
        a = x*(-4*M-4*m*M-4*l*M)
        b = x*(2*M+2*n*M+2*l*M+2*m*M)
        c = x*(-6-10*m-8*m*l-3*M-10*l-6*l*l-5*m*m*M-3*m*M-6*m*m
               +5*K*l*M+2*K*l*n*M+3*K*M+2*K*n*M+2*l*n*M+2*n*M
               -3*l*M+2*K*m*n*M+4*m*l*M+5*K*m*M+3*K*m*m*M+4*K*l*m*M
               -5*l*l*M+3*K*l*l*M+2*m*n*M)
    elif delta == (1,0,-2):
        c = n*(K*M-5*M-2)*(l+1)*(n-1)/(2*M)
    elif delta == (0,1,-2):
        c = n*(K*M-5*M-2)*(n-1)*(m+1)/(2*M)
    elif delta == (-2,0,1):
        c = l*(K*M-5*M-2)*(l-1)*(n+1)/M
    elif delta == (0,-2,1):
        c = m*(K*M-5*M-2)*(n+1)*(m-1)/M
    elif delta == (1,-2,0):
        c = 2*m*(K*M-M-2)*(l+1)*(m-1)/M
    elif delta == (-2,1,0):
        c = 2*l*(K*M-M-2)*(l-1)*(m+1)/M
    elif delta == (1,-1,-1):
        c = 2*m*n*(K*M-M-2)*(l+1)/M
    elif delta == (-1,1,-1):
        c = 2*l*n*(K*M-M-2)*(m+1)/M
    elif delta == (-1,-1,1):
        c = 2*m*l*(K*M-M-2)*(n+1)/M
    elif delta == (-2,0,0):
        x = -l*(l-1)/M
        a = x*(-4*M)
        c = x*(2*n*M+4*n+2*K*n*M+4*K*m*M+M+2+3*K*M)
    elif delta == (0,-2,0):
        x = -m*(m-1)/M
        a = x*(-4*M)
        c = x*(4*K*l*M+M+2*K*n*M+2*n*M+4*n+2*3*K*M)
    elif delta == (0,0,-2):
        x = -n*(n-1)/M
        b = x
        c = x*(K*m*M-m*M-l*M-2+K*M-2*l-2*m+K*l*M-M)
    elif delta == (-1,-1,0):
        x = -2*l*m/M
        a = x*(-4*M)
        b = x*(2*M)
        c = x*(K*M+2*K*n*M+4*K*m*M+4*K*l*M+4*n-M-2*n*M-2*l*M+2-2*m*M)
    elif delta == (-1,0,-1):
        x = -2*l*n/M
        a = x*(-2*M)
        b = x
        c = x*(K*M+2*K*l*M+2*K*m*M+K*n*M-4*m-2*l*M-2-4*l)
    elif delta == (0,-1,-1):
        x = -2*m*n/M
        a = x*(-2*M)
        b = x
        c = x*(2*K*l*M+K*n*M+K*M+2*K*m*M-2*m*M-4*l-2-4*m)
    elif delta == (-2,-1,0):
        c = 2*l*m*(K-1)*(l-1)
    elif delta == (-2,0,-1):
        c = l*n*(K*M-M-2)*(l-1)/M
    elif delta == (-1,-2,0):
        c = 2*m*l*(K-1)*(m-1)
    elif delta == (-1,0,-2):
        c = l*n*(K*M-M-2)*(n-1)/(2*M)
    elif delta == (0,-1,-2):
        c = m*n*(K*M-M-2)*(n-1)/(2*M)
    elif delta == (0,-2,-1):
        c = m*n*(K*M-M-2)*(m-1)/M
    elif delta == (-1,-1,-1):
        c = 2*l*m*n*(K*M-M-2)/M

    return (a,b,c)

def check_pekeris_H(Ho,So,H,S):
    ierr = 0
    for i in range(7):
        for j in range(i+1):
            if Ho[i,j] != H[i,j]:
                print("H mismatch ",i,j,Ho[i,j],H[i,j])
                ierr = ierr + 1
            if So[i,j] != S[i,j]:
                print("S mismatch ",i,j,So[i,j],S[i,j])
                ierr = ierr + 1
    print("There are a total of ",ierr," errors")
    return

def pekeris7x7(Z):
    # Tested and works for the 7x7 case. Will keep in program
    # for testing purposes.
    H = zeros((7,7),Float)
    S = zeros((7,7),Float)
    H[0,0] = -16.*Z + 5.
    S[0,0] = 16.
    H[1,0] = 4.*Z-4.
    S[1,0] = -4.
    H[2,0] = 28.*Z-6.
    S[2,0] = -28.
    H[3,0] = 1.
    S[3,0] = 0.
    H[4,0] = -4.*Z + 2.
    S[4,0] = 4.
    H[5,0] = -8.*Z
    S[5,0] = 8.
    H[6,0] = -4.*Z + 2.
    S[6,0] = 4.

    H[1,1] = -24.*Z + 15
    S[1,1] = 48.
    H[2,1] = -4.*Z+2.
    S[2,1] = -8.
    H[3,1] = 8.*Z-12.
    S[3,1] = -16.
    H[4,1] = 36.*Z - 10.
    S[4,1] = -60.
    H[5,1] = 0.
    S[5,1] = 8.

    H[2,2] = -112.*Z+26.
    S[2,2] = 144.
    H[3,2] = 0.
    S[3,2] = 4.
    H[4,2] = 16.*Z-12.
    S[4,2] = -16.
    H[5,2] = 88.*Z-12.
    S[5,2] = -104.
    H[6,2] = 44.*Z-14.
    S[6,2] = -72.

    H[3,3] = -32.*Z+31
    S[3,3] = 96.
    H[4,3] = -8.*Z+4.
    S[4,3] = -20.

    H[4,4] = -144.*Z + 54.
    S[4,4] = 336.
    H[5,4] = -8.*Z+4.
    S[5,4] = -32.
    H[6,4] = -4.*Z + 2.
    S[6,4] = -4.

    H[5,5] = -224.*Z + 34.
    S[5,5] = 320.
    H[6,5] = -16.*Z + 8.
    S[6,5] = 24.

    H[6,6] = -104.*Z + 25.
    S[6,6] = 208.

    H = symmetrize(H)
    S = symmetrize(S)

    return (H,S)

def symmetrize(M):
    rows,cols = M.shape
    if rows != cols:
        print( "Error: symmetrize only works for square matrices!")
        print( rows,cols)
        sys.exit()
    n = rows
    for i in range(n):
        for j in range(n):
            if j>i: M[i,j] = M[j,i]
    return M

def check_symmetry(M):
    rows,cols = M.shape
    if rows != cols: print( "check_symmetry: The matrix is not square!")
    for i in range(rows):
        for j in range(cols):
            if M[i,j] != M[j,i]:
                print( "check_symmetry: Mismatch",i,j,M[i,j],M[j,i])
    return

def testtransform():
    A = arange(9,typecode=Float)
    A = reshape(A,(3,3))
    B = ones((3,3),Float)
    print( transform(A,B))

def test_inv_sqrt():
    A = zeros((3,3),Float)
    A[0,0] = A[1,1] = A[2,2] = 3.
    A[1,2] = A[2,1] = A[0,1] = A[1,0] = 0.5
    print( A)
    X = inv_sqrt(A)
    print( transform(A,X))

