# PyQuante in Julia

# Experimenting with writing quantum chemistry in Julia
# This code is part of the PyQuante package, and is licensed under the same modified
# BSD license as the rest of PyQuante.

# Utility functions
factorial2(n) = prod(n:-2:1) # double factorial !!
dist2(dx,dy,dz) = dx*dx+dy*dy+dz*dz # Is there something in the standard library that does this?

function test_utilities()
    @assert factorial2(6)==48
end

# Basis function definitions
type PGBF
    expn::Float64
    x::Float64
    y::Float64
    z::Float64
    I::Int64
    J::Int64
    K::Int64
    norm::Float64
end

function pgbf(expn,x=0,y=0,z=0,I=0,J=0,K=0,norm=1)
    p = PGBF(expn,x,y,z,I,J,K,norm)
    normalize!(p)
    return p
end

function amplitude(bf::PGBF,x,y,z)
    dx,dy,dz = x-bf.x,y-bf.y,z-bf.z
    r2 = dist2(dx,dy,dz)
    return bf.norm*(dx^bf.I)*(dy^bf.J)*(dz^bf.K)*exp(-bf.expn*r2)
end

function normalize!(a::PGBF)
    olap = overlap(a,a)
    a.norm = 1/sqrt(olap)
end


type CGBF
    x::Float64
    y::Float64
    z::Float64
    I::Int64
    J::Int64
    K::Int64
    norm::Float64
    pgbfs::Array{PGBF,1}
    coefs::Array{Float64}
end

function cgbf(x=0,y=0,z=0,I=0,J=0,K=0)
    pgbfs = PGBF[]
    coefs = Int64[]
    c = CGBF(x,y,z,I,J,K,1.0,pgbfs,coefs)
    return c
end

function amplitude(bf::CGBF,x,y,z)
    s = 0
    for (c,pbf) in primitives(bf)
        s += c*amplitude(pbf,x,y,z)
    end
    return bf.norm*s
end

function normalize!(a::CGBF)
    olap = overlap(a,a)
    a.norm = 1/sqrt(olap)
end

primitives(a::CGBF) = zip(a.coefs,a.pgbfs)

function push!(cbf::CGBF,expn,coef)
    Base.push!(cbf.pgbfs,pgbf(expn,cbf.x,cbf.y,cbf.z,cbf.I,cbf.J,cbf.K))
    Base.push!(cbf.coefs,coef)
    normalize!(cbf)
end

function test_basis()
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @assert isapprox(amplitude(s,0,0,0),0.71270547)
    @assert isapprox(amplitude(px,0,0,0),0)
    c = cgbf(0.0,0.0,0.0)
    push!(c,1,1)
    @assert isapprox(amplitude(c,0,0,0),0.71270547)
end



# One-electron integrals
# Overlap matrix elements

function overlap(a::PGBF,b::PGBF)
    return a.norm*b.norm*overlap(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

function overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    gamma = aexpn+bexpn
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    rab2 = dist2(ax-bx,ay-by,az-bz) 
    pre = (pi/gamma)^1.5*exp(-aexpn*bexpn*rab2/gamma)
    wx = overlap1d(aI,bI,px-ax,px-bx,gamma)
    wy = overlap1d(aJ,bJ,py-ay,py-by,gamma)
    wz = overlap1d(aK,bK,pz-az,pz-bz,gamma)
    return pre*wx*wy*wz
end

function gaussian_product_center(a::PGBF,b::PGBF)
    return (a.expn*[a.x,a.y,a.z]+b.expn*[b.x,b.y,b.z])/(a.expn+b.expn)
end

function gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    return (aexpn*[ax,ay,az]+bexpn*[bx,by,bz])/(aexpn+bexpn)    
end

function overlap1d(la,lb,ax,bx,gamma)
    total = 0
    for i in 0:div(la+lb,2)
        total += binomial_prefactor(2i,la,lb,ax,bx)*factorial2(2i-1)/(2gamma)^i
    end
    return total
end

function binomial_prefactor(s,ia,ib,xpa,xpb)
    total = 0
    for t in 0:s
        if (s-ia) <= t <= ib
            total += binomial(ia,s-t)*binomial(ib,t)*xpa^(ia-s+t)*xpb^(ib-t)
        end
    end
    return total
end

function overlap(a::CGBF,b::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            s += ca*cb*overlap(abf,bbf)
        end
    end
    return a.norm*b.norm*s
end

function test_overlap()
    @assert overlap1d(0,0,0,0,1) == 1
    @assert gaussian_product_center(s,s) == [0,0,0]
    @assert isapprox(overlap(s,s),1)
    @assert isapprox(overlap(px,px),1)
    @assert isapprox(overlap(s,px),0)
    @assert binomial_prefactor(0,0,0,0,0) == 1
    @assert isapprox(overlap(c,c),1)
end


# Kinetic matrix elements
function kinetic(a::PGBF,b::PGBF)
    return a.norm*b.norm*kinetic(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

function kinetic(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    overlap0 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
    overlapx1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI+2,bJ,bK)
    overlapy1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ+2,bK)
    overlapz1 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK+2)
    overlapx2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-2,bJ,bK)
    overlapy2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-2,bK)
    overlapz2 = overlap(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-2)
    term0 = bexpn*(2*(bI+bJ+bK)+3)*overlap0
    term1 = -2*(bexpn^2)*(overlapx1+overlapy1+overlapz1)
    term2 = -0.5*(bI*(bI-1)*overlapx2+bJ*(bJ-1)*overlapy2+bK*(bK-1)*overlapz2)
    return term0+term1+term2
end

function kinetic(a::CGBF,b::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            s += ca*cb*kinetic(abf,bbf)
        end
    end
    return a.norm*b.norm*s
end

function test_kinetic()
    @assert isapprox(kinetic(1,0,0,0,0,0,0,1,0,0,0,0,0,0),2.9530518648229536)
    @assert isapprox(kinetic(s,s),1.5)
    @assert isapprox(kinetic(c,c),1.5)
end



# Nuclear attraction term
function Aterm(i,r,u,l1,l2,ax,bx,cx,gamma)
    term1 = (-1)^i*binomial_prefactor(i,l1,l2,ax,bx)
    term2 = (-1)^u*factorial(i)*cx^(i-2r-2u)
    term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
    return term1*term2*term3
end

function Aarray(l1,l2,a,b,c,g)
    Imax = l1+l2+1
    A = zeros(Float64,Imax)
    for i in 0:(Imax-1)
        for r in 0:div(i,2)
            for u in 0:div(i-2r,2)
                I = i-2r-u+1
                A[I] += Aterm(i,r,u,l1,l2,a,b,c,g)
            end
        end
    end
    return A
end

function nuclear_attraction(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,cx,cy,cz)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    gamma = aexpn+bexpn
    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcp2 = dist2(cx-px,cy-py,cz-pz)
    Ax = Aarray(aI,bI,px-ax,px-bx,px-cx,gamma)
    Ay = Aarray(aJ,bJ,py-ay,py-by,py-cy,gamma)
    Az = Aarray(aK,bK,pz-az,pz-bz,pz-cz,gamma)
    total = 0
    for I in 0:(aI+bI)
        for J in 0:(aJ+bJ)
            for K in 0:(aK+bK)
                total += Ax[I+1]*Ay[J+1]*Az[K+1]*Fgamma(I+J+K,rcp2*gamma)
            end
        end
    end
    return -2pi/gamma*exp(-aexpn*bexpn*rab2/gamma)*total
end

function nuclear_attraction(a::PGBF,b::PGBF,cx,cy,cz)
    return a.norm*b.norm*nuclear_attraction(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K,cx,cy,cz)
end

function Fgamma(m,x,SMALL=1e-12)
    #println("Fgamma($m,$x)")
    x = max(x,SMALL) # Evidently needs underflow protection
    return 0.5*x^(-m-0.5)*gammainc(m+0.5,x)
end

function gammainc(a,x)
    # This is the series version of gamma from pyquante. For reasons I don't get, it 
    # doesn't work around a=1. This works alright, but is only a stopgap solution
    # until Julia gets an incomplete gamma function programmed
    if abs(a-1) < 1e-3
        println("Warning: gammainc_series is known to have problems for a ~ 1")
    end
    gamser,gln = gser(a,x)
    return exp(gln)*gamser
end


function gser(a,x,ITMAX=100,EPS=3e-9)
    # Series representation of Gamma. NumRec sect 6.1.
    gln=lgamma(a)
    if x == 0
        return 0,gln
    end
    ap = a
    delt = s = 1/a
    for i in 1:ITMAX
        ap += 1
        delt *= (x/ap)
        s += delt
        if abs(delt) < abs(s)*EPS
            break
        end
    end
    return s*exp(-x+a*log(x)-gln),gln
end

function nuclear_attraction(a::CGBF,b::CGBF,cx,cy,cz)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            s += ca*cb*nuclear_attraction(abf,bbf,cx,cy,cz)
        end
    end
    return a.norm*b.norm*s
end


function test_nuke()
    @assert Aterm(0,0,0,0,0,0,0,0,0) == 1.0
    @assert Aarray(0,0,0,0,0,1) == [1.0]
    @assert Aarray(0,1,1,1,1,1) == [1.0, -1.0]
    @assert Aarray(1,1,1,1,1,1) == [1.5, -2.5, 1.0]
    @assert Aterm(0,0,0,0,0,0,0,0,1) == 1.0
    @assert Aterm(0,0,0,0,1,1,1,1,1) == 1.0
    @assert Aterm(1,0,0,0,1,1,1,1,1) == -1.0
    @assert Aterm(0,0,0,1,1,1,1,1,1) == 1.0
    @assert Aterm(1,0,0,1,1,1,1,1,1) == -2.0
    @assert Aterm(2,0,0,1,1,1,1,1,1) == 1.0
    @assert Aterm(2,0,1,1,1,1,1,1,1) == -0.5
    @assert Aterm(2,1,0,1,1,1,1,1,1) == 0.5

    # gammainc test functions. Test values taken from Mathematica
    # println("a=0.5 test")
    @assert maximum([gammainc(0.5,x) for x in 0:10]
                    -{0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968, 
                      1.77151, 1.77213, 1.77234, 1.77241, 1.77244}) < 1e-5
    @assert maximum([gammainc(1.5,x) for x in 0:10]
                    -{0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968, 
                      1.77151, 1.77213, 1.77234, 1.77241, 1.77244}) < 1e-5
    @assert maximum([gammainc(2.5,x) for x in 0:10]
                    -{0, 0.200538, 0.59898, 0.922271, 1.12165, 1.22933, 
                      1.2831, 1.30859, 1.32024, 1.32542, 1.32768}) < 1e-5

    @assert isapprox(nuclear_attraction(s,s,0,0,0),-1.59576912)
    @assert isapprox(nuclear_attraction(c,c,0,0,0),-1.595769)
end


# Two electron integrals

function coulomb(aexpn,ax,ay,az,aI,aJ,aK,
    bexpn,bx,by,bz,bI,bJ,bK,
    cexpn,cx,cy,cz,cI,cJ,cK,
    dexpn,dx,dy,dz,dI,dJ,dK)
    # This is the slow method of computing integrals from Huzinaga et al.
    # Use the HRR/VRR scheme from Head-Gordon & Pople instead

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,cx,cy,cz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)
    g1 = aexpn+bexpn
    g2 = cexpn+dexpn
    delta = 0.25*(1/g1+1/g2)
    
    Bx = Barray(aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)
    By = Barray(aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)
    Bz = Barray(aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)
    
    s = 0
    for I in 0:(aI+bI+cI+dI)
        for J in 0:(aJ+bJ+cJ+dJ)
            for K in 0:(aK+bK+cK+dK)
                s += Bx[I+1]*By[J+1]*Bz[K+1]*Fgamma(I+J+K,0.25*rpq2/delta)
            end
        end
    end
    return 2*pi^(2.5)/(g1*g2*sqrt(g1+g2))*exp(-aexpn*bexpn*rab2/g1)*exp(-cexpn*dexpn*rcd2/g2)*s
end

function coulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*coulomb(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
    c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
    d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end

fB(i,l1,l2,p,a,b,r,g) = binomial_prefactor(i,l1,l2,p-a,p-b)*B0(i,r,g)
B0(i,r,g) = fact_ratio2(i,r)*(4g)^(r-i)
fact_ratio2(a,b) = factorial(a,b)/factorial(a-2b)

function Bterm(i1,i2,r1,r2,u,l1,l2,l3,l4,Px,Ax,Bx,Qx,Cx,Dx,gamma1,gamma2,delta)
    # THO eq. 2.22
    return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)*(-1)^i2*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)
           *(-1)^u*fact_ratio2(i1+i2-2*(r1+r2),u)
           *(Qx-Px)^(i1+i2-2*(r1+r2)-2*u)/delta^(i1+i2-2*(r1+r2)-u)
end

function Barray(l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta)
    Imax = l1+l2+l3+l4+1
    B = zeros(Int64,Imax)
    for i1 in 0:(l1+l2)
        for i2 in 0:(l3+l4)
            for r1 in 0:div(i1,2)
                for r2 in 0:div(i2,2)
                    for u in 0:(div(i1+i2,2)-r1-r2)
                        I = i1+i2-2*(r1+r2)-u
                        B[I+1] += Bterm(i1,i2,r1,r2,u,l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta)
                    end
                end
            end
        end
    end
    return B
end

function coulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            for (cc,cbf) in primitives(c)
                for (cd,dbf) in primitives(d)
                    s += ca*cb*cc*cd*coulomb(abf,bbf,cbf,dbf)
                end
            end
        end
    end
    return a.norm*b.norm*c.norm*d.norm*s
end

function test_two()
    @assert isapprox(coulomb(1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0, 1, 0,0,0, 0,0,0),4.373355)
    @assert isapprox(coulomb(s,s,s,s),1.128379)
    @assert isapprox(coulomb(c,c,c,c),1.128379)
end

# Basis Set Data
sto3g = {
    # H
    [('S',
      [(3.4252509099999999, 0.15432897000000001),
       (0.62391373000000006, 0.53532813999999995),
       (0.16885539999999999, 0.44463454000000002)])],
    # He
    [('S',
      [(6.3624213899999997, 0.15432897000000001),
       (1.1589229999999999, 0.53532813999999995),
       (0.31364978999999998, 0.44463454000000002)])],
    # Li
    [('S',
      [(16.119575000000001, 0.15432897000000001),
       (2.9362007000000001, 0.53532813999999995),
       (0.79465050000000004, 0.44463454000000002)]),
     ('S',
      [(0.63628969999999996, -0.099967230000000004),
       (0.14786009999999999, 0.39951282999999999),
       (0.048088699999999998, 0.70011546999999996)]),
     ('P',
      [(0.63628969999999996, 0.15591627),
       (0.14786009999999999, 0.60768372000000004),
       (0.048088699999999998, 0.39195739000000002)])],
    # Be
    [('S',
      [(30.167871000000002, 0.15432897000000001),
       (5.4951153000000001, 0.53532813999999995),
       (1.4871927, 0.44463454000000002)]),
     ('S',
      [(1.3148331, -0.099967230000000004),
       (0.3055389, 0.39951282999999999),
       (0.099370700000000006, 0.70011546999999996)]),
     ('P',
      [(1.3148331, 0.15591627),
       (0.3055389, 0.60768372000000004),
       (0.099370700000000006, 0.39195739000000002)])],
    # B
    [('S',
      [(48.791113000000003, 0.15432897000000001),
       (8.8873622000000001, 0.53532813999999995),
       (2.4052669999999998, 0.44463454000000002)]),
     ('S',
      [(2.2369561, -0.099967230000000004),
       (0.51982050000000002, 0.39951282999999999),
       (0.16906180000000001, 0.70011546999999996)]),
     ('P',
      [(2.2369561, 0.15591627),
       (0.51982050000000002, 0.60768372000000004),
       (0.16906180000000001, 0.39195739000000002)])],
    # C
    [('S',
      [(71.616837000000004, 0.15432897000000001),
       (13.045095999999999, 0.53532813999999995),
       (3.5305122, 0.44463454000000002)]),
     ('S',
      [(2.9412493999999998, -0.099967230000000004),
       (0.68348310000000001, 0.39951282999999999),
       (0.22228990000000001, 0.70011546999999996)]),
     ('P',
      [(2.9412493999999998, 0.15591627),
       (0.68348310000000001, 0.60768372000000004),
       (0.22228990000000001, 0.39195739000000002)])],
    # N
    [('S',
      [(99.106168999999994, 0.15432897000000001),
       (18.052312000000001, 0.53532813999999995),
       (4.8856602000000002, 0.44463454000000002)]),
     ('S',
      [(3.7804559000000002, -0.099967230000000004),
       (0.87849659999999996, 0.39951282999999999),
       (0.28571439999999998, 0.70011546999999996)]),
     ('P',
      [(3.7804559000000002, 0.15591627),
       (0.87849659999999996, 0.60768372000000004),
       (0.28571439999999998, 0.39195739000000002)])],
    # O
    [('S',
      [(130.70931999999999, 0.15432897000000001),
       (23.808861, 0.53532813999999995),
       (6.4436083000000002, 0.44463454000000002)]),
     ('S',
      [(5.0331513000000001, -0.099967230000000004),
       (1.1695960999999999, 0.39951282999999999),
       (0.38038899999999998, 0.70011546999999996)]),
     ('P',
      [(5.0331513000000001, 0.15591627),
       (1.1695960999999999, 0.60768372000000004),
       (0.38038899999999998, 0.39195739000000002)])],
    # F
    [('S',
      [(166.67912999999999, 0.15432897000000001),
       (30.360811999999999, 0.53532813999999995),
       (8.2168206999999995, 0.44463454000000002)]),
     ('S',
      [(6.4648032000000004, -0.099967230000000004),
       (1.5022812000000001, 0.39951282999999999),
       (0.48858849999999998, 0.70011546999999996)]),
     ('P',
      [(6.4648032000000004, 0.15591627),
       (1.5022812000000001, 0.60768372000000004),
       (0.48858849999999998, 0.39195739000000002)])],
    # Ne
    [('S',
       [(207.01561000000001, 0.15432897000000001),
        (37.708151000000001, 0.53532813999999995),
        (10.205297, 0.44463454000000002)]),
      ('S',
       [(8.2463151000000003, -0.099967230000000004),
        (1.9162661999999999, 0.39951282999999999),
        (0.62322929999999999, 0.70011546999999996)]),
      ('P',
       [(8.2463151000000003, 0.15591627),
        (1.9162661999999999, 0.60768372000000004),
            (0.62322929999999999, 0.39195739000000002)])]
}

basis_set_data = {"sto3g" => sto3g}

# Atoms and Molecules
type Atom
    atno::Int64
    x::Float64
    y::Float64
    z::Float64
end

type Molecule
    atomlist::Array{Atom,1}
end

function push!(mol::Molecule,at::Atom)
    Base.push!(atomlist,at)
end

tobohr(x::Float64) = x/0.52918
function tobohr!(at::Atom)
    at.x /= 0.52918
    at.y /= 0.52918
    at.z /= 0.52918
end
function tobohr!(mol::Molecule)
    for at in mol.atomlist
        tobohr!(at)
    end
end

# Other molecule methods to implement
# nat, nuclear_repulsion, nel, nocc, nclosed, nopen, nup, ndown, stoich, mass,
# center_of_mass, center!

# Array of symbols, masses

# Sample molecules for tests
h2 = Molecule([Atom(1,  0.00000000,     0.00000000,     0.36628549),
               Atom(1,  0.00000000,     0.00000000,    -0.36628549)])

h2o = Molecule([Atom(8,   0.00000000,     0.00000000,     0.04851804),
                Atom(1,   0.75300223,     0.00000000,    -0.51923377),
                Atom(1,  -0.75300223,     0.00000000,    -0.51923377)])

ch4 = Molecule([Atom(6,   0.00000000,     0.00000000,     0.00000000),
                Atom(1,   0.62558332,    -0.62558332,     0.62558332),
                Atom(1,  -0.62558332,     0.62558332,     0.62558332),
                Atom(1,   0.62558332,     0.62558332,    -0.62558332),
                Atom(1,  -0.62558332,    -0.62558332,    -0.62558332)])

c6h6 = Molecule([ Atom(6,  0.98735329,     0.98735329,     0.00000000),
                  Atom(6,  1.34874967,    -0.36139639,     0.00000000),
                  Atom(6,  0.36139639,    -1.34874967,     0.00000000),
                  Atom(6, -0.98735329,    -0.98735329,     0.00000000),
                  Atom(6, -1.34874967,     0.36139639,     0.00000000),
                  Atom(6, -0.36139639,     1.34874967,     0.00000000),
                  Atom(1,  1.75551741,     1.75551741,     0.00000000),
                  Atom(1,  2.39808138,    -0.64256397,     0.00000000),
                  Atom(1,  0.64256397,    -2.39808138,     0.00000000),
                  Atom(1, -1.75551741,    -1.75551741,     0.00000000),
                  Atom(1, -2.39808138,     0.64256397,     0.00000000),
                  Atom(1, -0.64256397,     2.39808138,     0.00000000)])

# Convert to atomic units (bohr)
tobohr!(h2)
tobohr!(h2o)
tobohr!(ch4)
tobohr!(c6h6)

type BasisSet # list of CGBFs
    bfs::Array{CGBF,1}
end

basisset() = BasisSet(CGBF[])

function push!(basis::BasisSet,cbf::CGBF)
    Base.push!(basis.bfs,cbf)
end

function build_basis(mol::Molecule,name="sto3g")
    data = basis_set_data[name]
    basis_set = basisset()
    for atom in mol.atomlist
        for btuple in data[atom.atno]
            println(btuple)
            sym,primlist = btuple
            for (I,J,K) in sym2power[sym]
                cbf = cgbf(at.x,at.y,at.z,I,J,K)
                push!(basis_set,cbf)
                for (expn,coef) in primlist
                    push!(cbf,expn,coef)
                end
            end
        end
    end
    return basis_set
end

sym2power = {
    'S' => [(0,0,0)],
    'P' => [(1,0,0),(0,1,0),(0,0,1)],
    'D' => [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(1,0,1),(0,1,1)]
    } 

# Test functions
function test()
    test_utilities()
    test_basis()
    test_overlap()
    test_kinetic()
    test_nuke()
    test_two()
end

