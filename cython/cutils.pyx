def hello():
    return "Hello, World!"


def isnear(a,b,tol=1e-9): return abs(a-b)<tol

def fact2(n):
    """
    fact2(n) - n!!, double factorial of n
    >>> fact2(0)
    1
    >>> fact2(1)
    1
    >>> fact2(3)
    3
    >>> fact2(8)
    384
    >>> fact2(-1)
    1
    """
    return reduce(int.__mul__,xrange(n,0,-2),1)

