# This is a simplified version of the PyQuante1 and pyquante2 dft loops to work out
#  the bugs when implementing pyquante2.

import numpy as np

def pyq1_rohf(atomtuples=[(2,(0,0,0))],basis = '6-31G**',maxit=10,mult=3):
    from PyQuante import Ints,settings,Molecule
    from PyQuante.hartree_fock import get_energy
    from PyQuante.MG2 import MG2 as MolecularGrid
    from PyQuante.LA2 import mkdens,geigh,trace2,simx
    from PyQuante.Ints import getJ,getK
    
    print ("PyQ1 ROHF run")
    atoms = Molecule('Pyq1',atomlist=atomtuples,multiplicity=mult)

    bfs = Ints.getbasis(atoms,basis=basis)
    S,h,I2e = Ints.getints(bfs,atoms)

    nbf = norbs = len(bfs)
    nel = atoms.get_nel()

    nalpha,nbeta = atoms.get_alphabeta()

    enuke = atoms.get_enuke()
    orbe,orbs = geigh(h,S)
    eold = 0

    for i in range(maxit):
        Da = mkdens(orbs,0,nalpha)
        Db = mkdens(orbs,0,nbeta)
        Ja = getJ(I2e,Da)
        Jb = getJ(I2e,Db)
        Ka = getK(I2e,Da)
        Kb = getK(I2e,Db)

        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        eone = (trace2(Da,h) + trace2(Db,h))/2
        etwo = (trace2(Da,Fa) + trace2(Db,Fb))/2
        energy = (energya+energyb)/2 + enuke
        print (i,energy,eone,etwo,enuke)
        if abs(energy-eold) < 1e-5: break
        eold = energy

        Fa = simx(Fa,orbs)
        Fb = simx(Fb,orbs)
        # Building the approximate Fock matrices in the MO basis
        F = 0.5*(Fa+Fb)
        K = Fb-Fa

        # The Fock matrix now looks like
        #      F-K    |  F + K/2  |    F
        #   ---------------------------------
        #    F + K/2  |     F     |  F - K/2
        #   ---------------------------------
        #       F     |  F - K/2  |  F + K

        # Make explicit slice objects to simplify this
        do = slice(0,nbeta)
        so = slice(nbeta,nalpha)
        uo = slice(nalpha,norbs)
        F[do,do] -= K[do,do]
        F[uo,uo] += K[uo,uo]
        F[do,so] += 0.5*K[do,so]
        F[so,do] += 0.5*K[so,do]
        F[so,uo] -= 0.5*K[so,uo]
        F[uo,so] -= 0.5*K[uo,so]

        orbe,mo_orbs = np.linalg.eigh(F)
        orbs = np.dot(orbs,mo_orbs)
    return energy,orbe,orbs
def pyq1_rohf(atomtuples=[(2,(0,0,0))],basis = '6-31G**',maxit=10,mult=3):
    from PyQuante import Ints,settings,Molecule
    from PyQuante.hartree_fock import get_energy
    from PyQuante.MG2 import MG2 as MolecularGrid
    from PyQuante.LA2 import mkdens,geigh,trace2,simx
    from PyQuante.Ints import getJ,getK
    
    print ("PyQ1 ROHF run")
    atoms = Molecule('Pyq1',atomlist=atomtuples,multiplicity=mult)

    bfs = Ints.getbasis(atoms,basis=basis)
    S,h,I2e = Ints.getints(bfs,atoms)

    nbf = norbs = len(bfs)
    nel = atoms.get_nel()

    nalpha,nbeta = atoms.get_alphabeta()

    enuke = atoms.get_enuke()
    orbe,orbs = geigh(h,S)
    eold = 0

    for i in range(maxit):
        Da = mkdens(orbs,0,nalpha)
        Db = mkdens(orbs,0,nbeta)
        Ja = getJ(I2e,Da)
        Jb = getJ(I2e,Db)
        Ka = getK(I2e,Da)
        Kb = getK(I2e,Db)

        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        eone = (trace2(Da,h) + trace2(Db,h))/2
        etwo = (trace2(Da,Fa) + trace2(Db,Fb))/2
        energy = (energya+energyb)/2 + enuke
        print (i,energy,eone,etwo,enuke)
        if abs(energy-eold) < 1e-5: break
        eold = energy

        Fa = simx(Fa,orbs)
        Fb = simx(Fb,orbs)
        # Building the approximate Fock matrices in the MO basis
        F = 0.5*(Fa+Fb)
        K = Fb-Fa

        # The Fock matrix now looks like
        #      F-K    |  F + K/2  |    F
        #   ---------------------------------
        #    F + K/2  |     F     |  F - K/2
        #   ---------------------------------
        #       F     |  F - K/2  |  F + K

        # Make explicit slice objects to simplify this
        do = slice(0,nbeta)
        so = slice(nbeta,nalpha)
        uo = slice(nalpha,norbs)
        F[do,do] -= K[do,do]
        F[uo,uo] += K[uo,uo]
        F[do,so] += 0.5*K[do,so]
        F[so,do] += 0.5*K[so,do]
        F[so,uo] -= 0.5*K[so,uo]
        F[uo,so] -= 0.5*K[uo,so]

        orbe,mo_orbs = np.linalg.eigh(F)
        orbs = np.dot(orbs,mo_orbs)
    return energy,orbe,orbs
    

def pyq1_dft(atomtuples=[(2,(0,0,0))],basis = '6-31G**',maxit=10,
             xcname='SVWN'):
    from PyQuante import Ints,settings,Molecule
    from PyQuante.dft import getXC
    from PyQuante.MG2 import MG2 as MolecularGrid
    from PyQuante.LA2 import mkdens,geigh,trace2
    from PyQuante.Ints import getJ
    
    print ("PyQ1 DFT run")
    atoms = Molecule('Pyq1',atomlist=atomtuples)

    bfs = Ints.getbasis(atoms,basis=basis)
    S,h,Ints = Ints.getints(bfs,atoms)

    nclosed,nopen = nel//2,nel%2
    assert nopen==0
    enuke = atoms.get_enuke()

    grid_nrad = settings.DFTGridRadii
    grid_fineness = settings.DFTGridFineness

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness) 
    gr.set_bf_amps(bfs)

    orbe,orbs = geigh(h,S)
    eold = 0

    for i in range(maxit):
        D = mkdens(orbs,0,nclosed)
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,Vxc = getXC(gr,nel,functional=xcname)

        F = h+2*J+Vxc
        orbe,orbs = geigh(F,S)
        
        Ej = 2*trace2(D,J)
        Eone = 2*trace2(D,h)
        energy = Eone + Ej + Exc + enuke
        
        print (i,energy,Eone,Ej,Exc,enuke)
        if np.isclose(energy,eold):
            break
        eold = energy
    return energy

def func_compare():
    import matplotlib.pyplot as plt
    import pyquante2 as pyq2
    from PyQuante.DFunctionals import cvwn
    ns = np.linspace(0.,100)
    c2 = pyq2.dft.functionals.cvwn5(ns,ns)
    fc1 = [cvwn(n,n)[0] for n in ns]
    dfc1 = [cvwn(n,n)[1] for n in ns]
    plt.plot(ns,c2[0],label='f_vwn5/pyq2')#,marker='o',linestyle='None')
    plt.plot(ns,fc1,label='f_vwn5/pyq1')
    plt.plot(ns,c2[1],label='df_vwn5/pyq2')#,marker='o',linestyle='None')
    plt.plot(ns,dfc1,label='df_vwn5/pyq1')
    plt.show()

def pyq2_rohf(atomtuples=[(2,0,0,0)],basis = '6-31G**',maxit=10,xcname='svwn',
              mult=3):
    import pyquante2 as pyq2
    print ("pyq2 ROHF run")
    geo = pyq2.molecule(atomtuples,multiplicity=mult)
    bfs = pyq2.basisset(geo,name=basis)
    i1 = pyq2.onee_integrals(bfs,geo)
    i2 = pyq2.twoe_integrals(bfs)
    h = i1.T + i1.V
    orbe,orbs = pyq2.geigh(h,i1.S)
    eold = 0
    E0 = geo.nuclear_repulsion()
    nalpha,nbeta = geo.nup(),geo.ndown()
    norbs = len(bfs)

    for i in range(maxit):
        Da = pyq2.dmat(orbs,nalpha)
        Db = pyq2.dmat(orbs,nbeta)
        E1 = 0.5*pyq2.trace2(Da+Db,h)
        Ja,Ka = i2.get_j(Da),i2.get_k(Da)
        Jb,Kb = i2.get_j(Db),i2.get_k(Db)
        Fa = h + Ja + Jb - Ka
        Fb = h + Ja + Jb - Kb
        E2 = 0.5*(pyq2.trace2(Fa,Da)+pyq2.trace2(Fb,Db))
        energy = E0+E1+E2
        print (energy,E1,E2,E0)

        Fa = pyq2.utils.simx(Fa,orbs)
        Fb = pyq2.utils.simx(Fb,orbs)

        F = 0.5*(Fa+Fb)
        K = Fb-Fa
        # Make explicit slice objects to simplify this
        do = slice(0,nbeta)
        so = slice(nbeta,nalpha)
        uo = slice(nalpha,norbs)
        F[do,do] -= K[do,do]
        F[uo,uo] += K[uo,uo]
        F[do,so] += 0.5*K[do,so]
        F[so,do] += 0.5*K[so,do]
        F[so,uo] -= 0.5*K[so,uo]
        F[uo,so] -= 0.5*K[uo,so]

        E,cmo = np.linalg.eigh(F)
        orbs = np.dot(orbs,cmo)

    return

def pyq2_dft(atomtuples=[(2,0,0,0)],basis = '6-31G**',maxit=10,xcname='svwn'):
    import pyquante2 as pyq2
    print ("pyq2 DFT run")
    geo = pyq2.molecule(atomtuples)
    bfs = pyq2.basisset(geo,name=basis)
    i1 = pyq2.onee_integrals(bfs,geo)
    i2 = pyq2.twoe_integrals(bfs)
    grid = pyq2.grid(geo)
    h = i1.T + i1.V
    orbe,orbs = pyq2.geigh(h,i1.S)
    eold = 0
    grid.setbfamps(bfs)
    E0 = geo.nuclear_repulsion()

    for i in range(maxit):
        D = pyq2.dmat(orbs,geo.nocc())
        E1 = 2*pyq2.trace2(h,D)

        J = i2.get_j(D)
        Ej = 2*pyq2.trace2(J,D)

        Exc,Vxc = pyq2.get_xc(grid,0.5*D,xcname=xcname)

        energy = E0+E1+Ej+Exc
        F = h+2*J+Vxc

        orbe,orbs = pyq2.geigh(F,i1.S)

        print (i,energy,E1,Ej,Exc,E0)
        if np.isclose(energy,eold):
            break
        eold = energy
    return energy
    

if __name__ == '__main__':
    pyq1_rohf()
    pyq2_rohf()
    #pyq1_dft()
    #pyq2_dft()
    #func_compare()
    
