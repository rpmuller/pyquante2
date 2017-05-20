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

"""
