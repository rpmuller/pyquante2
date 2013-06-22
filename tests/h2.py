from pyquante2 import basisset,rhf,h2

def h2_orbs():
    bfs = basisset(h2,'sto3g')
    solver = rhf(h2,bfs)
    ens = solver.converge()
    print solver.orbs
    return

if __name__ == '__main__': h2_orbs()

