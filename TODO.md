# PyQuante 2 Todo List

## Soon
* Get He triplet XS dft test working. Probably involves changing the
iterator to the USCF one.
* Get LDA working

## Easy stuff to add:
* Get cvwn working
* Get clyp working
* Text cpbe
* Units for molecule._repr_html_
* __repr__ and _repr_html_ for basisset
* use numpy in shapes/viewer stuff
* from pyglet.gl import * ==> from pyglet import gl, gl.???
* get pyglet zooming to work, either through scrolling of modified clicks or something.
* Is there a redundancy in the Becke reweighting? Seems like I loop over all the atoms
  in the Ps loop, but also loop over all the atoms in becke_atomic_...

## Bug!!
I'm doing something dumb in my cute code that expands the cgbfs into
pgbfs in the integral code.  Everything should work fine if I'm doing
all pgbf or all cgbfs, but there's a bug if I have a mixture.

**Update**: Fixed for 1e code, not for 2e code.

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
* Should gaussian product center be in utils? should binomial whatever
* More intelligent way of handling physical constants, using one of 
  the many units packages that work with numpy (e.g. quantities, etc.)
* Make it possible for __call__ in basis functions to be a ufunc, so
  that we can evaluate wave functions along a set of points by just
  passing in a matrix??
* Should we require scipy to gain access to the incomplete gamma
  functions? Is there any way we can use the math.erf function in
  py27?
* How about scipy for the Legendre or Lebedev stuff?

## Things I'd like more of
* Itertools
* NamedTuples
* Cython
* Speed
* _repr_html_
* nosetests
* einsum
