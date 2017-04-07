The pyquante two-electron integral code uses three different methods:

* 'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa, S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966. http://doi.org/10.1143/JPSJ.21.2313 
* Augspurger, Bernholdt, and Dykstra. 'Concise, Open-Ended Implementation of Rys Polynomial Evaluation of Two-Electron Integrals.' J. Comp. Chem. 11 (8), 972-977 (1990). http://dx.doi.org/10.1002/jcc.540110809 
* M Head-Gordon and J Pople, ‘A method for two‐electron Gaussian integral and integral derivative evaluation using recurrence relations’. J. Chen. Phys. 89, 5777 (1988). http://dx.doi.org/10.1063/1.455553

By default, the Head-Gordon Pople method is used, since that’s both the simplest and the fastest, but the other methods are available. All three methods are coded in both pure-python and C.
