# PyQuante 2 Todo List

## Bug!!
I'm doing something dumb in my cute code that expands the cgbfs into pgbfs in the integral code.
Everything should work fine if I'm doing all pgbf or all cgbfs, but there's a bug if I have a mixture.

To fix this, I think I need 'else' statements:

    if b.contracted:
        return sum(cb*S(pb,a) for cb,pb in b)
    else:
        return S(b,a)

I don't think this is a big deal for the 1e code, but I think it's a problem for the 2e code.

## Next steps to take:
* Settings to ConfigParser
* Wrap libxc

## Questions to answer
* Should gaussian product center be in utils? should binomial whatever
* Make it possible for __call__ in basis functions to be a ufunc, so that we can evaluate 
  wave functions along a set of points by just passing in a matrix??
* Should we require scipy to gain access to the incomplete gamma functions? Is there any way we
  can use the math.erf function in py27?
* Should we consider doing everything in Cython?

## Things I'd like more of
* Itertools
* NamedTuples
* Cython
* Speed


