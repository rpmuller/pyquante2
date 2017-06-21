import numpy as np
from pyquante2 import *
import matplotlib.pyplot as plt

#
# Compute RHF on N2 molecule for interatomic
# separation R in a range with N points.
#

N = 50
R_vec = np.linspace(0.5, 3.0, N)
E_vec = np.zeros((N,))
for k in range(N):
    R = R_vec[k]
    #print "Solving for R = %g ..." % (R)
    n2 = molecule([(7,0,0,-R/2),(7,0,0,R/2)],units='Angstrom')
    bfs = basisset(n2,'sto3g')
    solver = rhf(n2,bfs)
    ens = solver.converge()

    E_vec[k] = solver.energy

    
plt.figure()
plt.plot(R_vec, E_vec)
plt.show()




                                        
