# PyQuante 2 Todo List

## Easy stuff to add:
* Get cvwn working
* Get clyp working
* Text cpbe
* Units for molecule._repr_html_
* __repr__ and _repr_html_ for basisset
* use numpy in shapes/viewer stuff
* Is there a redundancy in the Becke reweighting? Seems like I loop over all the atoms
  in the Ps loop, but also loop over all the atoms in becke_atomic_...

## DFT errors
* The slater exchange has a fairly high error (2e-6) in the test suite. 
  Is there a bug here? The value occurs at high densities. Does not 
  seem to be an incorrect factor value.
* cvwn has a similar bug with the dfb value. The erroneous value
  occurs at 0 values of nb.

## Next steps to take:
* Wrap libxc
* Settings to ConfigParser

## Questions to consider
* Move xyz readers into IO module?
* Should we require scipy to gain access to the incomplete gamma
  functions? 
* How about scipy for the Legendre or Lebedev stuff?

## Things I'd like more of
* Itertools
* NamedTuples
* Cython
* Speed
* _repr_html_
* nosetests
* einsum
