
bfs = basisset(h2,'sto3g')
solver = rhf(h2,bfs)
ens = solver.converge()
