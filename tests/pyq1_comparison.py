# This is a simplified version of the PyQuante1 and pyquante2 dft loops to work out
#  the bugs when implementing pyquante2.

import numpy as np
import matplotlib.pyplot as plt

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

    nel = atoms.get_nel()
    nclosed,nopen = nel//2,nel%2
    assert nopen==0
    enuke = atoms.get_enuke()

    grid_nrad = settings.DFTGridRadii
    grid_fineness = settings.DFTGridFineness

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness) 
    gr.set_bf_amps(bfs)

    orbe,orbs = geigh(h,S)
    eold = 0

    print "grid weights"
    print gr.weights()
    
    for i in range(maxit):
        D = mkdens(orbs,0,nclosed)
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,Vxc = getXC(gr,nel,functional=xcname)
        print Exc
        print Vxc
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

def pyq2_dft(atomtuples=[(2,0,0,0)],basis = '6-31G**',maxit=10,xcname='svwn'):
    print ("pyq2 DFT run")
    import pyquante2 as pyq2
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

    print "grid weights"
    print grid.points[:,3]
    
    for i in range(maxit):
        D = pyq2.dmat(orbs,geo.nocc())
        E1 = 2*pyq2.trace2(h,D)

        print "nel on grid:",np.dot(grid.points[:,3],grid.getdens(D))

        J = i2.get_j(D)
        Ej = 2*pyq2.trace2(J,D)

        Exc,Vxc = pyq2.get_xc(grid,0.5*D,xcname=xcname)

        print Exc
        print Vxc

        energy = E0+E1+Ej+Exc
        F = h+2*J+Vxc

        orbe,orbs = pyq2.geigh(F,i1.S)

        print (i,energy,E1,Ej,Exc,E0)
        if np.isclose(energy,eold):
            break
        eold = energy
    return energy
        
    

if __name__ == '__main__':
    pyq1_dft()
    pyq2_dft()
    #func_compare()
    
