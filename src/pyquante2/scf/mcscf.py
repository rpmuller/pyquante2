#!/usr/bin/env python
"""\
mcscf.py - multiconfigurational self-consistent field methods

The methods defined in this module seek to optimize a set of orbitals
with the energy expression:

  E = 2 f_i h_ii + a_ij J_ij + b_ij K_ij

Here
  f_i   Occupations of orbital i
  a_ij  Coulombic coefficient for orbitals i,j
  b_ij  Exchange coefficient for orbitals i,j

  h_ij  One electron Hamiltonian
        Sum of kinetic t_ij and nuclear attraction v_ij terms
  J_ij  Coulomb interaction between orbitals i,j
  K_ij  Exchange interaction between orbitals i,j

  i,j   indices of orbitals
  m,n   indices of basis functions

With different choices for the f,a,b terms, this energy expression can
describe
  - Closed shell RHF
  - Open shell ROHF
  - Generalized valence bond MC-SCF
  - More general MC-SCFs, CAS-SCF, etc.
We will start with the first three of these, and generalize to the fourth
when all is working with those.

A good reference to the general method is:
The Self-Consistent Equations for Generalized Valence Bond and Open-Shell
Hartree-Fock Wave Functions
Frank W. Bobrowicz and William A. Goddard, III, in
Methods of Electronic Structure Theory, H.F. Schaeffer, III, ed.
Plenum Publishing Corp., 1977
http://www.wag.caltech.edu/publications/sup/pdf/108.pdf

For RHF, ROHF, and GVB wave functions, the f,a,b terms take a fairly
simple form:

  f_i = 1      If the orbital is closed-shell (doubly occupied)
  f_i = 1/2    If the orbital is open-shell (singly-occupied, high spin)
  f_i = c_i^2  If the orbital is part of a GVB pair with coefficient c_i

  a_ij = 2 f_i f_j
  b_ij = -f_i f_j

except that:
  b_ij = -1/2               If i and j are both open-shell
  a_ii = f_i, b_ii = 0      If i is a pair orbital
  a_ij = 0, b_ij = -c_i c_j If i and j are in the same pair


The module is currently written to only go as complex as a GVB
wave function. However, you can implement a general MC-SCF wave
function by implementing a more complicated routine to update
the CI coefficients. 
"""
import logging
import numpy as np
from pyquante2 import molecule
from pyquante2.geo.samples import h,h2,lih,li
from pyquante2.basis.basisset import basisset
from pyquante2.utils import geigh,ao2mo
from pyquante2.ints.integrals import onee_integrals, twoe_integrals

logger = logging.getLogger(__name__)

try:
    np.set_printoptions(legacy='1.13')
except TypeError:
    pass

def gvb(geo,npair=0,basisname='sto3g',maxiter=25,verbose=False,
        return_orbs=False, input_orbs=None):
    """\
    This is a trivial test for the gvb module, because other
    pyquante modules are simpler if you're doing closed shell rhf,
    and should give the same results.

    # -0.46658184546856041 from uhf/sto3g
    >>> gvb(h)     # doctest: +ELLIPSIS
    -0.4665818...

    #  -1.117099582955609 from rhf/sto3g
    >>> gvb(h2)    # doctest: +ELLIPSIS
    -1.117099...

    >>> gvb(lih,maxiter=5)   # doctest: +ELLIPSIS
    -7.86073...

    >>> gvb(li,maxiter=5)    # doctest: +ELLIPSIS
    -7.31552...

    >>> gvb(h2,npair=1)      # doctest: +ELLIPSIS
    -1.13730...
    """
    # Get the basis set and the integrals
    bfs = basisset(geo,basisname)
    i1 = onee_integrals(bfs,geo)
    i2 = twoe_integrals(bfs)
    h = i1.T + i1.V

    # Get a guess for the orbitals
    if input_orbs is not None:
        U = input_orbs
    else:
        E,U = geigh(h,i1.S)

    # Set the parameters based on the molecule
    nopen = geo.nopen()
    ncore = geo.nclosed() - npair
    nocc = ncore + nopen + 2*npair
    norb = len(bfs)
    virt = range(nocc,norb)
    orbs_per_shell = get_orbs_per_shell(ncore,nopen,npair)
    nsh = len(orbs_per_shell)
    shell = orbital_to_shell_mapping(ncore,nopen,npair)
    Enuke = geo.nuclear_repulsion()

    f,a,b = fab(ncore,nopen,npair)

    if verbose:
        np.set_printoptions(precision=4)
        logger.info("**** PyQuante GVB ****")
        logger.info(f"{geo}")
        logger.info("Nuclear repulsion energy: %.3f" % Enuke)
        logger.info("Basis set: %s" % basisname)
        logger.info("  ncore/open/pair: %d,%d,%d" % (ncore,nopen,npair))
        logger.info("  occ/bf/orb: %d,%d,%d" % (nocc,len(bfs),norb))
        for i in range(nsh):
            logger.info("Shell %d" % i)
            logger.info("  occupation = %.2f" % f[i])
            logger.info("  orbitals in shell %s" % orbs_per_shell[i])
            logger.info("  couplings to other shells %s" % zip(a[i,:],b[i,:]))
        logger.info(f"Starting guess at orbitals:\n{U}")
        logger.info("Shell array: %s" % shell)
        logger.info("****")

    Eold = 0
    for it in range(maxiter):
        # Make all of the density matrices:
        Ds = [dmat_gen(U,orbs) for orbs in orbs_per_shell]
        # Compute the required Hamiltonian matrices:
        Js = [i2.get_j(D) for D in Ds]
        Ks = [i2.get_k(D) for D in Ds]

        # Perform the ROTION step and compute the energy
        Eel,Eone,Uocc = ROTION(U[:,:nocc],h,Js,Ks,f,a,b,nocc,shell,
                               verbose=verbose)
        if nsh > 1:
            U[:,:nocc] = Uocc

        # Perform the OCBSE step
        U = OCBSE(U,h,Js,Ks,f,a,b,orbs_per_shell,virt)

        #E = Enuke+Eone+Etwo
        E = Enuke+Eel
        Etwo = Eel-2*Eone

        # Update CI coefs
        coeffs = update_gvb_ci_coeffs(Uocc,h,Js,Ks,f,a,b,ncore,nopen,npair,
                                      orbs_per_shell,verbose)
        f,a,b = fab(ncore,nopen,npair,coeffs)

        if verbose:
            logger.info("---- %d :  %10.4f %10.4f %10.4f %10.4f" % ((it+1),E,Enuke,Eone,Etwo))
        if np.isclose(E,Eold):
            if verbose:
                logger.info("Energy converged")
            break
        Eold = E
    else:
        logger.warning("Maximum iterations (%d) reached without convergence" % (maxiter))
    if return_orbs:
        return E,U
    return E

def recompute_energy(Uocc,h,Js,Ks,f,a,b,nocc,shell):
    """\
    This is a helper routine to compute the GVB/ROHF energy expression.
    This routine should not be used in production code, since these
    terms are computed in ROTION.
    Eel,Eone = recompute_energy(Uocc,h,Js,Ks,f,a,b,nocc,shell)
    """
    nsh = len(f)
    hmo = ao2mo(h,Uocc)
    Jmo = [ao2mo(J,Uocc) for J in Js]
    Kmo = [ao2mo(K,Uocc) for K in Ks]
    Eone = sum(f[shell[i]]*hmo[i,i] for i in range(nocc))
    Fmo = [f[i]*hmo + sum(a[i,j]*Jmo[j] + b[i,j]*Kmo[j] for j in range(nsh))
           for i in range(nsh)]
    Eel = Eone + sum(Fmo[shell[i]][i,i] for i in range(nocc))
    return Eel,Eone

def ROTION(Uocc,h,Js,Ks,f,a,b,nocc,shell,verbose=False):
    """\
    Eel,Eone,Uocc = ROTION(Uocc,h,Js,Ks,f,a,b,nocc,shell,verbose)
    """
    nsh = len(f)
    hmo = ao2mo(h,Uocc)
    Jmo = [ao2mo(J,Uocc) for J in Js]
    Kmo = [ao2mo(K,Uocc) for K in Ks]
    Eone = sum(f[shell[i]]*hmo[i,i] for i in range(nocc))
    Fmo = [f[i]*hmo + sum(a[i,j]*Jmo[j] + b[i,j]*Kmo[j] for j in range(nsh))
           for i in range(nsh)]
    Eel = Eone + sum(Fmo[shell[i]][i,i] for i in range(nocc))
        
    Delta = np.zeros((nocc,nocc),'d')
    for i in range(nocc):
        ish = shell[i]
        for j in range(i):
            jsh = shell[j]
            if ish == jsh: continue
            # ish is now guaranteed to be larger than 0
            Jij = Jmo[ish][j,j]
            Kij = Kmo[ish][j,j]
            Gij = 2*(a[ish,ish]+a[jsh,jsh]-2*a[ish,jsh])*Kij \
                  + (b[ish,ish]+b[jsh,jsh]-2*b[ish,jsh])*(Jij+Kij)

            D0 = -(Fmo[jsh][i,j]-Fmo[ish][i,j])/\
                 (Fmo[jsh][i,i]-Fmo[ish][i,i]-Fmo[jsh][j,j]+Fmo[ish][j,j]\
                  +Gij)
            Delta[i,j] = D0
            Delta[j,i] = -D0
    if verbose:
        logger.info("ROTION Delta Matrix")
        logger.info(f"{Delta}")
    if nsh > 1:
        eD = expm(Delta)
        Uocc = np.dot(Uocc,eD)
    return Eel,Eone,Uocc

def OCBSE(U,h,Js,Ks,f,a,b,orbs_per_shell,virt):
    """\
    U = OCBSE(U,h,Js,Ks,f,a,b,orbs_per_shell,virt)

    Perform an Orthogonality Constrained Basis Set Expansion
    to mix the occupied orbitals with the virtual orbitals.
    See Bobrowicz/Goddard Sect 5.1.
    """
    Unew = np.zeros(U.shape,'d')
    nsh = len(f)
    for i,orbs in enumerate(orbs_per_shell):
        space = list(orbs) + list(virt)
        Fi = f[i]*h + sum(a[i,j]*Js[j]+b[i,j]*Ks[j] for j in range(nsh))
        Fi = ao2mo(Fi,U[:,space])
        Ei,Ci = np.linalg.eigh(Fi)
        Ui = np.dot(U[:,space],Ci)

        Unew[:,space] = Ui
    return Unew

def expm(M,tol=1e-6,maxit=30):
    """\
    Good time to cite 'Nineteen Dubious Ways to Compute the
    Exponential of a Matrix', Moler and Van Loan.
    http://epubs.siam.org/doi/abs/10.1137/S00361445024180
    >>> expm(np.zeros((2,2),'d'))
    array([[ 1.,  0.],
           [ 0.,  1.]])
    >>> expm(0.1*np.ones((2,2),'d'))
    array([[ 1.11070138,  0.11070138],
           [ 0.11070138,  1.11070138]])
    """
    n,m = M.shape
    assert n==m
    factor = 1.0
    X =  np.identity(n,'d')
    eM = np.identity(n,'d')
    for j in range(1,maxit):
        factor /= j
        X = np.dot(X,M)
        if np.linalg.norm(X) < tol:
            break
        eM += factor*X
    else:
        logger.warning("expm remainder = \n%s" % X)
        raise Exception("Maximum iterations reached in expm")
    return eM

def orbital_to_shell_mapping(ncore,nopen,npair):
    """\
    Map the orbitals to shells. All the core orbitals are in
    the first shell. Then each orbital has its own shell.
    >>> orbital_to_shell_mapping(1,0,0)
    [0]
    >>> orbital_to_shell_mapping(2,0,0)
    [0, 0]
    >>> orbital_to_shell_mapping(2,1,0)
    [0, 0, 1]
    """
    shell = [0 for i in range(ncore)]
    ncoreshells = 1 if ncore else 0
    for i in range(nopen+2*npair):
        shell.append(i+ncoreshells)
    return shell

def dmat_gen(c,indices): return np.dot(c[:,indices],c[:,indices].T)

def get_orbs_per_shell(ncore,nopen,npair):
    nocc = ncore+nopen+2*npair
    orbs_per_shell = []
    if ncore:
        orbs_per_shell.append(range(ncore))
    for i in range(ncore,nocc):
        orbs_per_shell.append([i])
    return orbs_per_shell

def update_gvb_ci_coeffs(Uocc,h,Js,Ks,f,a,b,ncore,nopen,npair,orbs_per_shell,
                         verbose=False):
    """\
    coeffs = update_gvb_ci_coeffs(Uocc,h,Js,Ks,f,a,b,ncore,nopen,npair,
                                  orbs_per_shell,verbose=False)
    """
    # consider reusing the transformed MO integral elements from ROTION
    #  (or moving the module there)
    coeffs = np.zeros((2*npair,),'d')
    nsh = len(f)
    ncoresh = 1 if ncore else 0
    hmo = ao2mo(h,Uocc)
    Jmo = [ao2mo(J,Uocc) for J in Js]
    Kmo = [ao2mo(K,Uocc) for K in Ks]

    for i in range(npair):
        ish = ncoresh+nopen+i
        jsh = ish+1
        iorb = orbs_per_shell[ish][0]
        jorb = orbs_per_shell[jsh][0]
        H11 = 2*hmo[iorb,iorb] + Jmo[ish][iorb,iorb]
        H22 = 2*hmo[jorb,jorb] + Jmo[jsh][jorb,jorb]
        K12 = Kmo[ish][jorb,jorb] # == Kmo[jsh][iorb,iorb] (checked)
        for k in range(nsh):
            if k == ish or k == jsh: continue
            H11 += f[k]*(2*Jmo[ksh][iorb,iorb]-Kmo[ksh][iorb,iorb])
            H22 += f[k]*(2*Jmo[ksh][jorb,jorb]-Kmo[ksh][jorb,jorb])

        H = np.zeros((2,2),'d')
        H[0,1] = H[1,0] = K12
        H[0,0] = H11
        H[1,1] = H22
        if verbose:
            logger.info("GVB CI Matrix for pair %d\n%s" % (i,H))
        E,C = np.linalg.eigh(H)
        if verbose:
            logger.info("GVB CI Eigenvector for pair %d\n%s" % (i,C))
            logger.info("GVB CI Eigenvalues for pair %d\n%s" % (i,E))
        coeffs[2*i:(2*i+2)] = C[0]
    return coeffs

def guess_gvb_ci_coeffs(npair,state='0'):
    """
    Make a guess at the CI coefficients for the GVB pairs.
    The orbitals are ordered:
      (pair 1, natural orbital 1),
      (pair 2, natural orbital 1),
      ...
      (pair 1, natural orbital 2),
      (pair 2, natural orbital 2),
      ...
    whereas the gvb ci coefficients are ordered:
      p1n1,p1n2,p2n1,p2n2,...
    >>> guess_gvb_ci_coeffs(1)
    array([ 1.,  0.])
    >>> guess_gvb_ci_coeffs(2)
    array([ 1.,  0.,  1.,  0.])
    """
    inv_rt2 = 1/np.sqrt(2)
    coeffs = []
    for i in range(npair):
        if state == '1':
            coeffs.append(0)
            coeffs.append(1)
        elif state == '+':
            coeffs.append(inv_rt2)
            coeffs.append(inv_rt2)
        elif state == '-':
            coeffs.append(inv_rt2)
            coeffs.append(-inv_rt2)            
        else: # Default state=0
            coeffs.append(1)
            coeffs.append(0)
    return np.array(coeffs,'d')

def fab(ncore,nopen,npair,coeffs=None):
    """\
    Create arrays over shells for the coefficients of the
    energy expression:
      E = 2 f_i h_ii + a_ij J_ij + b_ij K_ij
    Here
      f_i   Occupations of orbital i
      a_ij  Coulombic coefficient for orbitals i,j
      b_ij  Exchange coefficient for orbitals i,j

    Shells in GVB are a little confusing. All core orbitals are
    in a single shell, and then each occupied orbital has its own
    shell.

    Closed-shell, e.g. H2
    >>> f,a,b = fab(1,0,0)
    >>> f
    array([ 1.])
    >>> a
    array([[ 2.]])
    >>> b
    array([[-1.]])

    Single open-shell orbital, e.g. H
    >>> f,a,b = fab(0,1,0)
    >>> f
    array([ 0.5])
    >>> a
    array([[ 0.]])
    >>> b
    array([[ 0.]])

    Single gvb pair, e.g. H2
    >>> f,a,b = fab(0,0,1)
    >>> f
    array([ 1.,  0.])
    >>> a
    array([[ 1.,  0.],
           [ 0.,  0.]])
    >>> b
    array([[ 0., -0.],
           [-0.,  0.]])

    Single gvb pair, e.g. H2, explicitly specifying gvb coefficients
    >>> rt12 = np.sqrt(0.5)
    >>> f,a,b = fab(0,0,1,[rt12,rt12])
    >>> f
    array([ 0.5,  0.5])
    >>> a
    array([[ 0.5,  0. ],
           [ 0. ,  0.5]])
    >>> b
    array([[ 0. , -0.5],
           [-0.5,  0. ]])

    This tests a special case for a single open shell, where a[i,i] = 
    b[i,i] = 0 for the open shell
    >>> f,a,b = fab(1,1,0)
    >>> f
    array([ 1. ,  0.5])
    >>> a
    array([[ 2.,  1.],
           [ 1.,  0.]])
    >>> b
    array([[-1. , -0.5],
           [-0.5,  0. ]])
    """
    ncoresh = 1 if ncore else 0
    nsh = ncoresh+nopen+2*npair
    f = np.zeros(nsh,'d')
    a = np.zeros((nsh,nsh),'d')
    b = np.zeros((nsh,nsh),'d')

    if npair > 0 and coeffs is None:
        coeffs = guess_gvb_ci_coeffs(npair,'0')

    # f array
    if ncore:
        f[0] = 1
    for i in range(nopen):
        f[ncoresh+i] = 0.5

    # This is a little tricky:
    # Assume that the coeffs are arranged p1n1,p1n2,p2n1,p2n2,...
    # But the orbitals are arranged by occupation, p1n1,p2n1,...,p1n2,p2n2,...
    # I may rethink this -- seems unnaturally complex
    for p in range(npair):
        i = ncoresh+nopen+p
        f[i] = coeffs[2*p]**2
        f[i+npair] = coeffs[2*p+1]**2

    # Basic a,b assumptions:
    for i in range(nsh):
        for j in range(nsh):
            a[i,j] = 2*f[i]*f[j]
            b[i,j] = -f[i]*f[j]
            
    # Corrections
    # b_ij = -1/2               If i and j are both open-shell
    for i in range(ncoresh,ncoresh+nopen):
        for j in range(ncoresh,ncoresh+nopen):
            b[i,j] = -0.5

    # If there is only one open shell i, a[i,i]=b[i,i] = 0
    if nopen == 1:
        iop = ncoresh+nopen-1
        a[iop,iop] = b[iop,iop] = 0
    # a_ii = f_i, b_ii = 0      If i is a pair orbital
    for i in range(ncoresh+nopen,ncoresh+nopen+2*npair):
        a[i,i] = f[i]
        b[i,i] = 0
    # a_ij = 0, b_ij = -c_i c_j If i and j are in the same pair
    for i in range(ncoresh+nopen,ncoresh+nopen+npair):
        a[i,i+npair] = a[i+npair,i] = 0
        b[i,i+npair] = b[i+npair,i] = -coeffs[i]*coeffs[i+npair]
    return f,a,b

if __name__ == '__main__':
    #import doctest; doctest.testmod()
    #gvb(h,maxiter=5,verbose=True)   # -0.46658
    #gvb(lih,maxiter=5,verbose=True)   # -7.86073
    #gvb(li,maxiter=5,verbose=True)   # -7.3155
    #gvb(h2,npair=1,verbose=True) # RHF = -1.1171, GVB ?= -1.1373
    #
    # h- example for Rajib: Does not currently work with GVB.
    # RHF energy is -0.4224, GVB energy is -0.3964
    h_m = molecule(atomlist = [(1,0,0,0)],charge=-1,name="H-")
    from pyquante2.orbman import orbman
    E,orbs = gvb(h_m,maxiter=10,basisname='6-31g',npair=0,verbose=True,
               return_orbs=True)
    o_hf = orbman(orbs,basisset(h_m,'6-31g'),h_m)
    E,orbs = gvb(h_m,maxiter=10,basisname='6-31g',npair=1,verbose=True,
               return_orbs=True,input_orbs=orbs)
    o_gvb = orbman(orbs,basisset(h_m,'6-31g'),h_m)

    for i in [0,1]:
        logger.info("Orbital %d after HF" % i)
        o_hf[i]
        logger.info("  GVB")
        o_gvb[i]
    
