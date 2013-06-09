# PyQuante 2 Todo List

## Next steps to take:
* Allow contr_hrr to be called by fixing the
    Cannot convert Python object to 'double *'
  error.
* Settings to ConfigParser
* HF example
* Wrap libxc

## Questions to answer
* Should gaussian product center be in utils? should binomial whatever
* Make it possible for __call__ in basis functions to be a ufunc, so that we can evaluate 
  wave functions along a set of points by just passing in a matrix??
* Should we require scipy to gain access to the incomplete gamma functions? Is there any way we
  can use the math.erf function in py27?
* Should we consider doing everything in Cython?

