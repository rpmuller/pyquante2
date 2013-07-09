cimport cints

def S(a,b):
    if b.contracted:
        return sum(cb*S(pb,a) for (cb,pb) in b)
    elif a.contracted:
        return sum(ca*S(b,pa) for (ca,pa) in a)
    return a.norm*b.norm*cints.overlap(a.exponent,a.powers[0],a.powers[1],a.powers[2],
                                       a.origin[0],a.origin[1],a.origin[2],
                                       b.exponent,b.powers[0],b.powers[1],b.powers[2],
                                       b.origin[0],b.origin[1],b.origin[2])

def T(a,b):
    if b.contracted:
        return sum(cb*T(pb,a) for (cb,pb) in b)
    elif a.contracted:
        return sum(ca*T(b,pa) for (ca,pa) in a)
    return a.norm*b.norm*cints.kinetic(a.exponent,a.powers[0],a.powers[1],a.powers[2],
                                       a.origin[0],a.origin[1],a.origin[2],
                                       b.exponent,b.powers[0],b.powers[1],b.powers[2],
                                       b.origin[0],b.origin[1],b.origin[2])

def V(a,b,C):
    if b.contracted:
        return sum(cb*V(pb,a,C) for (cb,pb) in b)
    elif a.contracted:
        return sum(ca*V(b,pa,C) for (ca,pa) in a)
    return cints.nuclear_attraction(
        a.origin[0],a.origin[1],a.origin[2],a.norm,a.powers[0],
        a.powers[1],a.powers[2],a.exponent,
        b.origin[0],b.origin[1],b.origin[2],b.norm,b.powers[0],
        b.powers[1],b.powers[2],b.exponent,
        C[0],C[1],C[2])
