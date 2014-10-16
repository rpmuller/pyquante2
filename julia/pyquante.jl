# coding: utf-8
"""
# PyQuante in Julia
Experimenting with writing quantum chemistry in Julia

## Utility functions
"""

factorial2(n) = prod(n:-2:1) # double factorial !!
dist2(dx,dy,dz) = dx*dx+dy*dy+dz*dz # Is there something in the standard library that does this?

function pairs(n::Int64,which="diag")
    function _it()
        for j in 1:n
            start = j+1
            if which=="diag"
                start = j
            elseif which=="rect"
                start = 1
            end
            for i in start:n
                produce((j,i))
            end
        end
    end
    Task(_it)
end

triangle(i::Int64) = div(i*(i+1),2)
triangle(i::Int64,j::Int64) = i<j ? triangle(j-1)+i : triangle(i-1)+j

function iiterator(n::Int64)
    function _it()
        for (i,j) in pairs(n)
            ij = triangle(i-1,j)
            for (k,l) in pairs(n)
                kl = triangle(k-1,l)
                if ij <= kl
                    produce((i,j,k,l))
                end
            end
        end
    end
    Task(_it)
end

iindex(i::Int64,j::Int64,k::Int64,l::Int64) = triangle(triangle(i,j),triangle(k,l)) 

trace2(A,B) = sum(A.*B)

function test_utils()
    @assert factorial2(6)==48
    @assert collect(pairs(3)) == {(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)}
    @assert collect(pairs(3,"subdiag")) == {(1,2),(1,3),(2,3)}
    @assert collect(pairs(2,"rect")) == {(1,1),(1,2),(2,1),(2,2)}
    @assert iindex(1,1,1,1) == 1
    @assert iindex(1,1,1,2) == iindex(1,1,2,1) == iindex(1,2,1,1) == iindex(2,1,1,1) == 2
    @assert iindex(1,1,2,2) == iindex(2,2,1,1) == 4
end


# ## Basis function definitions

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

function normalize!(pbf::PGBF)
    pbf.norm /= sqrt(overlap(pbf,pbf))
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
    coefs::Array{Float64,1}
end

cgbf(x=0,y=0,z=0,I=0,J=0,K=0) = CGBF(x,y,z,I,J,K,1.0,PGBF[],Float64[])

function amplitude(bf::CGBF,x,y,z)
    s = 0
    for (c,pbf) in primitives(bf)
        s += c*amplitude(pbf,x,y,z)
    end
    return bf.norm*s
end

function normalize!(bf::CGBF)
    bf.norm /= sqrt(overlap(bf,bf))
end

primitives(a::CGBF) = zip(a.coefs,a.pgbfs)

function push!(cbf::CGBF,expn,coef)
    Base.push!(cbf.pgbfs,pgbf(expn,cbf.x,cbf.y,cbf.z,cbf.I,cbf.J,cbf.K))
    Base.push!(cbf.coefs,coef)
    normalize!(cbf)
end

function contract(f,a::CGBF,b::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            s += ca*cb*f(abf,bbf)
        end
    end
    return a.norm*b.norm*s
end

function contract(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
        for (cb,bbf) in primitives(b)
            for (cc,cbf) in primitives(c)
                for (cd,dbf) in primitives(d)
                    s += ca*cb*cc*cd*f(abf,bbf,cbf,dbf)
                end
            end
        end
    end
    return a.norm*b.norm*c.norm*d.norm*s
end


function test_pgbf()
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @assert isapprox(amplitude(s,0,0,0),0.71270547)
    @assert isapprox(amplitude(px,0,0,0),0)
end
function test_cgbf()
    c = cgbf(0.0,0.0,0.0)
    push!(c,1,1)
    @assert isapprox(amplitude(c,0,0,0),0.71270547)
    c2 = cgbf(0,0,0)
    push!(c2,1,0.2)
    push!(c2,0.5,0.2)
    @assert isapprox(overlap(c2,c2),1)
end

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

"""
## One-electron integrals
### Overlap matrix elements
"""

function overlap(a::PGBF,b::PGBF)
    return a.norm*b.norm*overlap(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
    b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

overlap(a::CGBF,b::CGBF) = contract(overlap,a,b)

function overlap(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                 aI::Int64,aJ::Int64,aK::Int64, 
                 bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
                 bI::Int64,bJ::Int64,bK::Int64)
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

function gaussian_product_center(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                                    bexpn::Float64,bx::Float64,by::Float64,bz::Float64)
    return (aexpn*[ax,ay,az]+bexpn*[bx,by,bz])/(aexpn+bexpn)    
end

function overlap1d(la::Int64,lb::Int64,ax::Float64,bx::Float64,gamma::Float64)
    total = 0
    for i in 0:div(la+lb,2)
        total += binomial_prefactor(2i,la,lb,ax,bx)*factorial2(2i-1)/(2gamma)^i
    end
    return total
end

function binomial_prefactor(s::Int64,ia::Int64,ib::Int64,xpa::Float64,xpb::Float64)
    #println("binomial_prefactor($s,$ia,$ib,$xpa,$xpb)")
    total = 0
    for t in 0:s
        if (s-ia) <= t <= ib
            total += binomial(ia,s-t)*binomial(ib,t)*xpa^(ia-s+t)*xpb^(ib-t)
        end
    end
    return total
end


function test_overlap()
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @assert overlap1d(0,0,0.,0.,1.) == 1
    @assert gaussian_product_center(s,s) == [0,0,0]
    @assert isapprox(overlap(s,s),1)
    @assert isapprox(overlap(px,px),1)
    @assert isapprox(overlap(s,px),0)
    @assert binomial_prefactor(0,0,0,0.,0.) == 1
end

"""
### Kinetic matrix elements
"""

function kinetic(a::PGBF,b::PGBF)
    return a.norm*b.norm*kinetic(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                                b.expn,b.x,b.y,b.z,b.I,b.J,b.K)
end

function kinetic(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
            bexpn::Float64,bx::Float64,by::Float64,bz::Float64,bI::Int64,bJ::Int64,bK::Int64)
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

kinetic(a::CGBF,b::CGBF) = contract(kinetic,a,b)


function test_kinetic()
    s = pgbf(1.0)
    c = cgbf(0.0,0.0,0.0)
    push!(c,1,1)
    @assert isapprox(amplitude(c,0,0,0),0.71270547)
    @assert isapprox(kinetic(1.,0.,0.,0.,0,0,0,1.,0.,0.,0.,0,0,0),2.9530518648229536)
    @assert isapprox(kinetic(s,s),1.5)
    @assert isapprox(kinetic(c,c),1.5)
end

"""
### Nuclear attraction term
"""
function Aterm(i::Int64,r::Int64,u::Int64,l1::Int64,l2::Int64,ax::Float64,bx::Float64,cx::Float64,gamma::Float64)
    term1 = (-1)^(i+u)*binomial_prefactor(i,l1,l2,ax,bx)
    term2 = factorial(i)*cx^(i-2r-2u)
    term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
    return term1*term2*term3
end

function Aarray(l1::Int64,l2::Int64,a::Float64,b::Float64,c::Float64,g::Float64)
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

function nuclear_attraction(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                            aI::Int64,aJ::Int64,aK::Int64,
                            bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
                            bI::Int64,bJ::Int64,bK::Int64,
                            cx::Float64,cy::Float64,cz::Float64)
    #print("na($aexpn,$ax,$ay,$az,$aI,$aJ,$aK,$bexpn,$bx,$by,$bz,$bI,$bJ,$bK,$cx,$cy,$cz)=")
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
    val=-2pi*exp(-aexpn*bexpn*rab2/gamma)*total/gamma
    #println(val)
    #println((Ax,Ay,Az,rcp2*gamma,Fgamma(0,rcp2*gamma)))
    return val
end

function nuclear_attraction(a::PGBF,b::PGBF,cx::Float64,cy::Float64,cz::Float64)
    return a.norm*b.norm*nuclear_attraction(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
                                            b.expn,b.x,b.y,b.z,b.I,b.J,b.K,cx,cy,cz)
end
nuclear_attraction(a::PGBF,b::PGBF,c::Atom) = c.atno*nuclear_attraction(a,b,c.x,c.y,c.z)
nuclear_attraction(a::PGBF,b::PGBF,m::Molecule) = sum([nuclear_attraction(a,b,c) for c in m.atomlist])

function Fgamma(m::Int64,x::Float64,SMALL::Float64=1e-12)
    #println("Fgamma($m,$x)")
    x = max(x,SMALL) # Evidently needs underflow protection
    return 0.5*x^(-m-0.5)*gammainc(m+0.5,x)
end

function gammainc(a::Float64,x::Float64)
    # This is the series version of gamma from pyquante. For reasons I don't get, it 
    # doesn't work around a=1. This works alright, but is only a stopgap solution
    # until Julia gets an incomplete gamma function programmed
    if abs(a-1) < 1e-3
        println("Warning: gammainc_series is known to have problems for a ~ 1")
    end
    if x < (a+1.0)
        #Use the series representation
        gam,gln = gser(a,x)
    else 
        #Use continued fractions
        gamc,gln = gcf(a,x)
        gam = 1-gamc
    end
    return exp(gln)*gam
end

function gser(a::Float64,x::Float64,ITMAX::Int64=100,EPS::Float64=3e-9)
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

function gcf(a::Float64,x::Float64,ITMAX::Int64=200,EPS::Float64=3e-9,FPMIN::Float64=1e-30)
    #Continued fraction representation of Gamma. NumRec sect 6.1"
    gln=lgamma(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    for i in 1:ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if abs(d) < FPMIN
            d=FPMIN
        end
        c=b+an/c
        if abs(c) < FPMIN
            c=FPMIN
        end
        d=1./d
        delt=d*c
        h=h*delt
        if abs(delt-1.) < EPS
            break
        end
    end
    gammcf = exp(-x+a*log(x)-gln)*h
    return gammcf,gln
end


# Need a nested scope to squeeze this into the contract function
function nuclear_attraction(a::CGBF,b::CGBF,cx::Float64,cy::Float64,cz::Float64)
    na(a,b) = nuclear_attraction(a,b,cx,cy,cz)
    contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,c::Atom)
    na(a,b) = nuclear_attraction(a,b,c)
    contract(na,a,b)
end
function nuclear_attraction(a::CGBF,b::CGBF,m::Molecule)
    na(a,b) = nuclear_attraction(a,b,m)
    contract(na,a,b)
end

function test_a_terms()
    @assert Aterm(0,0,0,0,0,0.,0.,0.,0.) == 1.0
    @assert Aarray(0,0,0.,0.,0.,1.) == [1.0]
    @assert Aarray(0,1,1.,1.,1.,1.) == [1.0, -1.0]
    @assert Aarray(1,1,1.,1.,1.,1.) == [1.5, -2.5, 1.0]
    @assert Aterm(0,0,0,0,0,0.,0.,0.,1.) == 1.0
    @assert Aterm(0,0,0,0,1,1.,1.,1.,1.) == 1.0
    @assert Aterm(1,0,0,0,1,1.,1.,1.,1.) == -1.0
    @assert Aterm(0,0,0,1,1,1.,1.,1.,1.) == 1.0
    @assert Aterm(1,0,0,1,1,1.,1.,1.,1.) == -2.0
    @assert Aterm(2,0,0,1,1,1.,1.,1.,1.) == 1.0
    @assert Aterm(2,0,1,1,1,1.,1.,1.,1.) == -0.5
    @assert Aterm(2,1,0,1,1,1.,1.,1.,1.) == 0.5
end

function test_gamma()
    # gammainc test functions. Test values taken from Mathematica
    # println("a=0.5 test")
    @assert maximum([gammainc(0.5,float(x)) for x in 0:10]
            -{0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968, 
                1.77151, 1.77213, 1.77234, 1.77241, 1.77244}) < 1e-5

    # println("a=1.5 test")
    @assert maximum([gammainc(1.5,float(x)) for x in 0:10]
            -{0, 1.49365, 1.69181, 1.7471, 1.76416, 1.76968, 
                1.77151, 1.77213, 1.77234, 1.77241, 1.77244}) < 1e-5
    # println("a=2.5 test")
    @assert maximum([gammainc(2.5,float(x)) for x in 0:10]
            -{0, 0.200538, 0.59898, 0.922271, 1.12165, 1.22933, 
                1.2831, 1.30859, 1.32024, 1.32542, 1.32768}) < 1e-5
end

function test_na()
    s = pgbf(1.0)
    c = cgbf(0.0,0.0,0.0)
    push!(c,1,1)
    @assert isapprox(amplitude(c,0,0,0),0.71270547)
    @assert isapprox(nuclear_attraction(s,s,0.,0.,0.),-1.59576912)
    @assert isapprox(nuclear_attraction(c,c,0.,0.,0.),-1.595769)
end

function test_fgamma()
    for (x,res) in {(0.,1),
                    (30.,0.161802),
                    (60.,0.1144114),
                    (90.,0.0934165),
                    (120.,0.08090108),
                    (300.,0.051166336)}
        @assert isapprox(res,Fgamma(0,x))
    end
end

#todo make into a test
function test_one()
    s1 = pgbf(1)
    s2 = pgbf(1,0,1,0)
    x=y=0.
    println("S: $(overlap(s1,s2))")
    println("T: $(kinetic(s1,s2))")
    for z in linspace(0,1,5)
        println("V: $z $(nuclear_attraction(s1,s2,x,y,z))")
    end
end


function test_na2()
    li,h = lih.atomlist
    bfs = build_basis(lih)
    s1,s2,x,y,z,h1s = bfs.bfs
    @assert isapprox(nuclear_attraction(s1,s1,lih),-8.307532656)
end

"""
## Two electron integrals
"""

function coulomb(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,
                 aI::Int64,aJ::Int64,aK::Int64,
                 bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
                 bI::Int64,bJ::Int64,bK::Int64,
                 cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,
                 cI::Int64,cJ::Int64,cK::Int64,
                 dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,
                 dI::Int64,dJ::Int64,dK::Int64)
    # This is the slow method of computing integrals from Huzinaga et al.
    # Use the HRR/VRR scheme from Head-Gordon & Pople instead

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)

    g1 = aexpn+bexpn
    g2 = cexpn+dexpn
    delta = 0.25*(1/g1+1/g2)
    
    Bx = Barray(aI,bI,cI,dI,px,ax,bx,qx,cx,dx,g1,g2,delta)
    By = Barray(aJ,bJ,cJ,dJ,py,ay,by,qy,cy,dy,g1,g2,delta)
    Bz = Barray(aK,bK,cK,dK,pz,az,bz,qz,cz,dz,g1,g2,delta)
    
    s = 0
    #println("$(aI+bI+cI+dI),$(aJ+bJ+cJ+dJ),$(aK+bK+cK+dK)")
    for I in 0:(aI+bI+cI+dI)
        for J in 0:(aJ+bJ+cJ+dJ)
            for K in 0:(aK+bK+cK+dK)
                #println("coul: $I,$J,$K,$(Bx[I+1]),$(By[J+1]),$(Bz[K+1])")
                s += Bx[I+1]*By[J+1]*Bz[K+1]*Fgamma(I+J+K,0.25*rpq2/delta)
            end
        end
    end
    return 2pi^(2.5)/(g1*g2*sqrt(g1+g2))*exp(-aexpn*bexpn*rab2/g1)*exp(-cexpn*dexpn*rcd2/g2)*s
end

function coulomb(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*coulomb(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
        b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
        c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
        d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end

fB(i::Int64,l1::Int64,l2::Int64,p::Float64,a::Float64,b::Float64,r::Int64,g::Float64) = binomial_prefactor(i,l1,l2,p-a,p-b)*B0(i,r,g)
B0(i::Int64,r::Int64,g::Float64) = fact_ratio2(i,r)*(4g)^(r-i)
fact_ratio2(a::Int64,b::Int64) = factorial(a,b)/factorial(a-2b)

function Bterm(i1::Int64,i2::Int64,r1::Int64,r2::Int64,u::Int64,
               l1::Int64,l2::Int64,l3::Int64,l4::Int64,
               Px::Float64,Ax::Float64,Bx::Float64,Qx::Float64,Cx::Float64,Dx::Float64,
               gamma1::Float64,gamma2::Float64,delta::Float64)
    # THO eq. 2.22
    #print("Bterm($i1,$i2,$r1,$r2,$u,$l1,$l2,$l3,$l4,$Px,$Ax,$Bx,$Qx,$Cx,$Dx,$gamma1,$gamma2,$delta)=")
    val = (-1)^(i2+u)*fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)*(
            fact_ratio2(i1+i2-2*(r1+r2),u)*(Qx-Px)^(i1+i2-2*(r1+r2)-2*u)/delta^(i1+i2-2*(r1+r2)-u))
    #println("$val")
    return val
end

function Barray(l1::Int64,l2::Int64,l3::Int64,l4::Int64,p::Float64,a::Float64,b::Float64,
                q::Float64,c::Float64,d::Float64,g1::Float64,g2::Float64,delta::Float64)
    Imax = l1+l2+l3+l4+1
    B = zeros(Float64,Imax)
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

coulomb(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb,a,b,c,d)


function test_two_terms()
    @assert fB(0,0,0,0.0,0.0,0.0,0,2.0) == 1
    @assert fB(0,0,0,1.0,1.0,1.0,0,2.0) == 1
    @assert fB(0,0,0,0.0,0.0,0.0,0,2.0 ) == 1
    @assert fB(1,0,1,0.0,0.0,0.0,0,2.0 ) == 0.125
    @assert B0(0,0,2.0) == 1
    @assert fact_ratio2(0,0) == 1
    @assert Bterm(0,0,0,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==1
    @assert Bterm(0,1,0,0,0,0,0,0,1,0.0,0.0,0.0,0.0,0.0,0.0,2.0,2.0,0.25)==0
end


function test_coul1()
    s = pgbf(1.0)
    px = pgbf(1.0,0,0,0,1,0,0)
    @assert coulomb(s,s,s,px)==0 # 0
    @assert isapprox(coulomb(s,s,px,px), 0.9403159725793305 )
end

function coulomb_hgp(a::PGBF,b::PGBF,c::PGBF,d::PGBF)
    return a.norm*b.norm*c.norm*d.norm*hrr(a.expn,a.x,a.y,a.z,a.I,a.J,a.K,
        b.expn,b.x,b.y,b.z,b.I,b.J,b.K,
        c.expn,c.x,c.y,c.z,c.I,c.J,c.K,
        d.expn,d.x,d.y,d.z,d.I,d.J,d.K)
end
coulomb_hgp(a::CGBF,b::CGBF,c::CGBF,d::CGBF) = contract(coulomb_hgp,a,b,c,d)

function hrr(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
    bexpn::Float64,bx::Float64,by::Float64,bz::Float64,bI::Int64,bJ::Int64,bK::Int64,
    cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
    dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,dI::Int64,dJ::Int64,dK::Int64)
    if bI > 0
        return hrr(aexpn,ax,ay,az,aI+1,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
        (ax-bx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI-1,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif bJ > 0
        return hrr(aexpn,ax,ay,az,aI,aJ+1,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
            (ay-by)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ-1,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif bK > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK+1,bexpn,bx,by,bz,bI,bJ,bK-1,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK) +
            (az-bz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK-1,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
    elseif dI > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI+1,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK) +
            (cx-dx)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI-1,dJ,dK)
    elseif dJ > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ+1,cK,dexpn,dx,dy,dz,dI,dJ-1,dK) +
            (cy-dy)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ-1,dK)
    elseif dK > 0
        return hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK+1,dexpn,dx,dy,dz,dI,dJ,dK-1) +
            (cz-dz)*hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK-1)
    end
    return vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,0)
end


function vrr(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
        bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
        cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
        dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,m::Int64)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
    #println("P: $px,$py,$pz, Q: $qx,$qy,$qz, W: $wx,$wy,$wz, $zeta,$eta")
    
    val = 0
    if cK>0
        val = (qz-cz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m) +
            (wz-qz)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val1=$val")
        if cK>1
            val += 0.5*(cK-1)/eta*(
                vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
                vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-2,dexpn,dx,dy,dz,m+1) )
        #println("val2=$val")
        end
        if aK>0
            val += 0.5*aK/(zeta+eta)*
                vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                    cexpn,cx,cy,cz,cI,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        end
        #println("val3=$val")
        return val
    elseif cJ>0
        val = (qy-cy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m) +
        (wy-qy)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-1,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val4=$val")
        if cJ>1
            val += 0.5*(cJ-1)/eta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-2,cK,dexpn,dx,dy,dz,m+1)
            )
        #println("val5=$val")
        end
        if aJ>0
            val += 0.5*aJ/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ-1,cK,dexpn,dx,dy,dz,m+1)
        end
        #println("val6=$val")
        return val
    elseif cI>0
        val = (qx-cx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) +
        (wx-qx)*vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK-1,dexpn,dx,dy,dz,m+1)
        #println("val7=$val")
        if cI>1
            val += 0.5*(cI-1)/eta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m) -
            zeta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-2,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val8=$val")
        if aI>0
            val += 0.5*aI/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
        end
        #println("val9=$val")
        return val
    elseif aK>0
        val = (pz-az)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
        (wz-pz)*vrr(aexpn,ax,ay,az,aI,aJ,aK-1,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val10=$val")
        if aK>1
            val += 0.5*(aK-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ,aK-2,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI-1,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val11=$val")
        return val
    elseif aJ>0
        val = (py-ay)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m)+
        (wy-py)*vrr(aexpn,ax,ay,az,aI,aJ-1,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val12=$val")
        if aJ>1
            val += 0.5*(aJ-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI,aJ-2,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val13=$val")
        return val
    elseif aI>0
        val = (px-ax)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) +
        (wx-px)*vrr(aexpn,ax,ay,az,aI-1,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
        #println("val14=$val")
        if aI>1
            val += 0.5*(aI-1)/zeta*(
            vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m) -
            eta/(zeta+eta)*
            vrr(aexpn,ax,ay,az,aI-2,aJ,aK,bexpn,bx,by,bz,
                cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,m+1)
            )
        end
        #println("val15=$val")
        return val
    end

    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    #println("rab2=$rab2,rcd2=$rcd2,rpq2=$rpq2,T=$T,Kab=$Kab,Kcd=$Kcd")
    return Kab*Kcd/sqrt(zeta+eta)*Fgamma(m,T)
end


function vrr_iter(aexpn::Float64,ax::Float64,ay::Float64,az::Float64,aI::Int64,aJ::Int64,aK::Int64,
        bexpn::Float64,bx::Float64,by::Float64,bz::Float64,
        cexpn::Float64,cx::Float64,cy::Float64,cz::Float64,cI::Int64,cJ::Int64,cK::Int64,
        dexpn::Float64,dx::Float64,dy::Float64,dz::Float64,M::Int64)
    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    qx,qy,qz = gaussian_product_center(cexpn,cx,cy,cz,dexpn,dx,dy,dz)
    zeta,eta = aexpn+bexpn,cexpn+dexpn
    wx,wy,wz = gaussian_product_center(zeta,px,py,pz,eta,qx,qy,qz)
    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcd2 = dist2(cx-dx,cy-dy,cz-dz)
    rpq2 = dist2(px-qx,py-qy,pz-qz)
    T = zeta*eta/(zeta+eta)*rpq2
    Kab = sqrt(2)pi^1.25/zeta*exp(-aexpn*bexpn*rab2/zeta)
    Kcd = sqrt(2)pi^1.25/eta*exp(-cexpn*dexpn*rcd2/eta)
    mtot = aI+aJ+aK+cI+cJ+cK+M

    vrr_terms = zeros(Float64,(aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,mtot+1))
    
    for m in 0:mtot
        vrr_terms[1,1,1, 1,1,1, m+1] = Fgamma(m,T)*Kab*Kcd/sqrt(zeta+eta)
    end
    
    for i in 0:(aI-1)
        for m in 0:(mtot-i-1)
            vrr_terms[i+2,1,1, 1,1,1, m+1] = (
                 (px-ax)*vrr_terms[i+1,1,1, 1,1,1, m+1] + 
                 (wx-px)*vrr_terms[i+1,1,1, 1,1,1, m+2])

            if i>0
                vrr_terms[i+2,1,1, 1,1,1, m+1] += i/2/zeta*(
                    vrr_terms[i,1,1, 1,1,1, m+1] -
                    eta/(zeta+eta)*vrr_terms[i,1,1, 1,1,1, m+2])
            end
        end
    end
#=
    for j in 0:(aJ-1)
        for i in 0:aI
            for m in 0:(mtot-i-j-1)
                println(("b",i,j,m))
                vrr_terms[i+1,j+2,1, 1,1,1, m+1] = (
                (py-ay)*vrr_terms[i+1,j+1,1, 1,1,1, m+1] +
                (wy-py)*vrr_terms[i+1,j+1,1, 1,1,1, m+2])
                if j>0
                    vrr_terms[i+1,j+2,1, 1,1,1, m+1] += j/2/zeta*(
                        vrr_terms[i+1,j,1, 1,1,1, m+1] -
                        eta/(zeta+eta)*vrr_terms[i+1,j,1, 1,1,1, m+2])
                end
            end
        end
    end

    for k in 0:(aK-1)
        for j in 0:aJ
            for i in 0:aI
                for m in 0:(mtot-i-j-k-1)
                    println(("c",i,j,k,m))
                    vrr_terms[i+1,j+1,k+2, 1,1,1, m+1] = (
                    (pz-az)*vrr_terms[i+1,j+1,k+1, 1,1,1, m+1] +
                    (wz-pz)*vrr_terms[i+1,j+1,k+1, 1,1,1, m+2])
                    if k>0
                        vrr_terms[i+1,j+1,k+2, 1,1,1, m+1] += k/2/zeta*(
                        vrr_terms[i+1,j+1,k, 1,1,1, m+1] -
                        eta/(zeta+eta)*vrr_terms[i+1,j+1,k, 1,1,1, m+2])
                    end
                end
            end
        end
    end

    for q in 0:(cI-1)
        for k in 0:aK
            for j in 0:aJ
                for i in 0:aI
                    for m in 0:(mtot-i-j-k-q-1)
                        println(("d",i,j,k,q,m))
                        vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] = (
                        (qx-cx)*vrr_terms[i+1,j+1,k+1, q+1,1,1, m+1] +
                        (wx-qx)*vrr_terms[i+1,j+1,k+1, q+1,1,1, m+2])
                        if q>0
                            vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] += q/2/eta*(
                            vrr_terms[i+1,j+1,k+1, q,1,1, m+1] -
                            eta/(zeta+eta)*vrr_terms[i+1,j+1,k+1, q,1,1, m+2])
                        end
                        if i>0
                            vrr_terms[i+1,j+1,k+1, q+2,1,1, m+1] += (
                            i/2/(zeta+eta)*vrr_terms[i,j+1,j+1, q+1,1,1, m+2])
                        end
                    end
                end
            end
        end
    end

    for r in 0:(cJ-1)
        for q in 0:cI
            for k in 0:aK
                for j in 0:aJ
                    for i in 0:aI
                        for m in 0:(mtot-i-j-k-q-r-1)
                            println(("e",i,j,k,q,r,m))
                            vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] = (
                            (qy-cy)*vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+1] +
                            (wy-qy)*vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+2])
                            if r>0
                                vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] += r/2/eta*(
                                vrr_terms[i+1,j+1,k+1, q+1,r+1,1, m+1] -
                                zeta/(zeta+eta)*vrr_terms[i+1,j+1,k+1, q+1,r,1, m+2])
                            end
                            if j>0
                                vrr_terms[i+1,j+1,k+1, q+1,r+2,1, m+1] += (
                                j/2/(zeta+eta)*vrr_terms[i+1,j,k+1, q+1,r+1,1, m+2])
                            end
                        end
                    end
                end
            end
        end
    end

    for s in 0:(cK-1)
        for r in 0:cJ
            for q in 0:cI
                for k in 0:aK
                    for j in 0:aJ
                        for i in 0:aI
                            for m in 0:(mtot-i-j-k-q-r-s-1)
                                println(("f",i,j,k,q,r,s,m))
                                vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2, m+1] = (
                                (qz-cz)*vrr_terms[i+1,j+1,k+1, q+1,r+1,s+1, m+1]+
                                (wz-qz)*vrr_terms[i+1,j+1,k+1, q+1,r+1,s+1, m+2])
                                if s>0
                                    vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2,m+1] += s/2/eta*(
                                    vrr_terms[i+1,j+1,k+1,q+1,r+1,s, m+1] -
                                    zeta/(zeta+eta)*vrr_terms[i+1,j+1,k+1,q+1,r+1,s,m+2])
                                end
                                if k>0
                                    vrr_terms[i+1,j+1,k+1, q+1,r+1,s+2,m+1] += (
                                    k/2/(zeta+eta)*vrr_terms[i+1,j+1,k, q+1,r+1,s+1,m+2])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("before return")
    =#
    @show vrr_terms[aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,M+1]
    vrr_terms[aI+1,aJ+1,aK+1,cI+1,cJ+1,cK+1,M+1]
end

function test_vrr()
    ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
    aexpn=bexpn=cexpn=dexpn=1.0
    aI=aJ=aK=0
    cI=cJ=cK=0
    M=0

    for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in {
            (0.,0.,0., 0,0,0, 0,0,0, 4.37335456733),
            (0.,0.,0., 1,0,0, 1,0,0, 0.182223107579),
            (0.,0.,0., 0,1,0, 0,1,0, 0.182223107579),
            (0.,0.,0., 0,0,1, 0,0,1, 0.182223107579),

            (0.,0.,0., 2,0,0, 2,0,0,  0.223223306785),
            (0.,0.,0., 0,2,0, 0,2,0,  0.223223306785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.223223306785),

            (1.,2.,3., 1,0,0, 1,0,0, -5.63387712455e-06),
            (1.,2.,3., 0,1,0, 0,1,0, -0.000116463120359),
            (1.,2.,3., 0,0,1, 0,0,1, -0.000301178525749),

            (1.,2.,3., 2,0,0, 2,0,0, 0.000225033081978),
            (1.,2.,3., 0,2,0, 0,2,0, 0.000610247078796),
            (1.,2.,3., 0,0,2, 0,0,2, 0.00134278307956),

            (0.,0.,0., 1,1,0, 1,1,0, 0.0136667330685),
            (0.,0.,0., 0,1,1, 0,1,1, 0.0136667330685),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0136667330685),

            (3.,2.,1., 1,1,0, 1,1,0, 5.97677147819e-05),
            (3.,2.,1., 0,1,1, 0,1,1, 1.57429039496e-06),
            (3.,2.,1., 1,0,1, 1,0,1, 4.00292836291e-06)
        }

        val1 = vrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
        val2 = vrr(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
            aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
        @assert isapprox(val1,val2)
        @assert isapprox(val1,result)
        val3 = vrr_iter(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,M)
        val4 = vrr_iter(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,
            aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,M)
        @show (val1,val2,val3,val4)
    end
end

function test_hrr()
    ax=ay=az=bx=by=bz=cx=cy=cz=dx=dy=dz=0.0
    aexpn=bexpn=cexpn=dexpn=1.0
    aI=aJ=aK=0
    bI,bJ,bK = 1,0,1
    cI=cJ=cK=0
    dI,dJ,dK = 1,0,1


    for (ax,ay,az, aI,aJ,aK, cI,cJ,cK, result) in {
            (0.,0.,0., 0,0,0, 0,0,0, 0.0136667330685),
            (0.,0.,0., 1,0,0, 1,0,0, 0.00821630976139),
            (0.,0.,0., 0,1,0, 0,1,0, 0.00122024402397),
            (0.,0.,0., 0,0,1, 0,0,1, 0.00821630976139),

            (0.,0.,0., 2,0,0, 2,0,0,   0.0039759617781),
            (0.,0.,0., 0,2,0, 0,2,0,   0.000599953311785),
            (0.,0.,0., 0,0,2, 0,0,2,  0.0039759617781),

            (1.,2.,3., 1,0,0, 1,0,0, -1.18513079462e-06),
            (1.,2.,3., 0,1,0, 0,1,0,  -4.66999128258e-06),
            (1.,2.,3., 0,0,1, 0,0,1, -3.47437366868e-05),

            (1.,2.,3., 2,0,0, 2,0,0, 2.81002247462e-06),
            (1.,2.,3., 0,2,0, 0,2,0, 7.09856891538e-06),
            (1.,2.,3., 0,0,2, 0,0,2, 3.62153023224e-05),

            (0.,0.,0., 1,1,0, 1,1,0, 0.000599953311785),
            (0.,0.,0., 0,1,1, 0,1,1, 0.000599953311785),
            (0.,0.,0., 1,0,1, 1,0,1, 0.0116431617287),

            (3.,2.,1., 1,1,0, 1,1,0, 7.37307761485e-06),
            (3.,2.,1., 0,1,1, 0,1,1, 2.53332441858e-07),
            (3.,2.,1., 1,0,1, 1,0,1, 2.4521155336e-06)
        }
        #println("hrr($aexpn,$ax,$ay,$az,$aI,$aJ,$aK,$bexpn,$bx,$by,$bz,$bI,$bJ,$bK,")
        #println("    $cexpn,$cx,$cy,$cz,$cI,$cJ,$cK,$dexpn,$dx,$dy,$dz,$dI,$dJ,$dK)")
        val1 = hrr(aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK,
            cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK)
        val2 = hrr(cexpn,cx,cy,cz,cI,cJ,cK,dexpn,dx,dy,dz,dI,dJ,dK,
            aexpn,ax,ay,az,aI,aJ,aK,bexpn,bx,by,bz,bI,bJ,bK)
        #@show val1,val2,result
        @assert isapprox(val1,val2)
        @assert isapprox(val1,result)
    end
end


# ## Basis Set Data
# Note use of curly braces here. Julia assumes that if you have square braces, you want
# things flattened as much as possible (to be as fast as possible, I guess). Curlys 
# preserve the list structure the way I would expect from Python

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
basis_set_data = {"sto3g" => sto3g};


# ## Atoms and Molecules

nuclear_repulsion(a::Atom,b::Atom)= a.atno*b.atno/sqrt(dist2(a.x-b.x,a.y-b.y,a.z-b.z))
function nuclear_repulsion(mol::Molecule)
    nr = 0
    for (i,j) in pairs(nat(mol),"subdiag")
        nr += nuclear_repulsion(mol.atomlist[i],mol.atomlist[j])
    end
    return nr
end

nel(mol::Molecule) = sum([at.atno for at in mol.atomlist])
nat(mol::Molecule) = length(mol.atomlist)

# Other molecule methods to implement
# nocc, nclosed, nopen, nup, ndown, stoich, mass,
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

lih = Molecule([Atom(3,    0.00000000,     0.00000000,    -0.53999756),
                Atom(1,    0.00000000,     0.00000000,     1.08999756)])

# Convert to atomic units (bohr)
tobohr!(h2)
tobohr!(h2o)
tobohr!(ch4)
tobohr!(c6h6)
tobohr!(lih)

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
            sym,primlist = btuple
            for (I,J,K) in sym2power[sym]
                cbf = cgbf(atom.x,atom.y,atom.z,I,J,K)
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


function test_geo_basis()
    @assert isapprox(nuclear_repulsion(h2),0.7223600367)
    @assert nel(h2) == 2
    @assert nel(h2o) == 10
    @assert length(sto3g)==10
    bfs = build_basis(h2)
    @assert length(bfs.bfs)==2
    l,r = bfs.bfs
    @assert isapprox(overlap(l,l),1)
    @assert isapprox(overlap(r,r),1)
    @assert isapprox(overlap(l,r),0.66473625)
    @assert isapprox(kinetic(l,l),0.76003188)
    @assert isapprox(kinetic(r,r),0.76003188)
    @assert isapprox(kinetic(l,r),0.24141861181119084)
    @assert isapprox(coulomb(l,l,l,l), 0.7746059439196398)
    @assert isapprox(coulomb(r,r,r,r), 0.7746059439196398)
    @assert isapprox(coulomb(l,l,r,r), 0.5727937653511646)
    @assert isapprox(coulomb(l,l,l,r), 0.4488373301593464)
    @assert isapprox(coulomb(l,r,l,r), 0.3025451156654606)
    bfs = build_basis(h2o)

    s1,s2,px,py,pz,hl,hr = bfs.bfs
    @assert isapprox(coulomb(s1,s2,hl,hr),0.03855344493645537)
    @assert isapprox(coulomb(s1,pz,hl,hr),-0.0027720110485359053)
    @assert isapprox(coulomb(s1,hl,pz,hr),-0.010049491284827426)
    @assert coulomb(s1,py,hl,hr)==0
    @assert coulomb(s1,hl,py,hr)==0
end

function all_1e_ints(bfs::BasisSet,mol::Molecule)
    n = length(bfs.bfs)
    S = Array(Float64,(n,n))
    T = Array(Float64,(n,n))
    V = Array(Float64,(n,n))
    for (i,j) in pairs(n)
        a,b = bfs.bfs[i],bfs.bfs[j]
        S[i,j] = S[j,i] = overlap(a,b)
        T[i,j] = T[j,i] = kinetic(a,b)
        V[i,j] = V[j,i] = nuclear_attraction(a,b,mol)
    end
    return S,T,V
end

function all_twoe_ints(bflist,ERI=coulomb)
    n = length(bflist.bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = Array(Float64,totlen)
    for (i,j,k,l) in iiterator(n)
        ints2e[iindex(i,j,k,l)] = ERI(bflist.bfs[i],bflist.bfs[j],bflist.bfs[k],bflist.bfs[l])
    end
    return ints2e
end

function make2JmK(D::Array{Float64,2},Ints::Array{Float64,1})
    n = size(D,1)
    G = Array(Float64,(n,n))
    D1 = reshape(D,n*n)
    temp = Array(Float64,n*n)
    for (i,j) in pairs(n)
        kl = 1
        for (k,l) in pairs(n,"rect")
            temp[kl] = 2*Ints[iindex(i,j,k,l)]-Ints[iindex(i,k,j,l)]
            kl += 1
        end
        G[i,j] = G[j,i] = dot(D1,temp)
    end
    return G
end

dmat(U::Array{Float64,2},nocc::Int64) = U[:,1:nocc]*U[:,1:nocc]'


function rhf(mol::Molecule,MaxIter::Int64=8,verbose::Bool=false)
    bfs = build_basis(mol)
    S,T,V = all_1e_ints(bfs,mol)
    Ints = all_twoe_ints(bfs)
    h = T+V
    E,U = eig(h,S)
    Enuke = nuclear_repulsion(mol)
    nclosed,nopen = divrem(nel(mol),2)
    Eold = 0
    Energy = 0
    println("Nel=$(nel(mol)) Nclosed=$nclosed")
    if verbose
        println("S=\n$S")
        println("h=\n$h")
        println("T=\n$T")
        println("V=\n$V")
        println("E: $E")
        println("U: $U")
        println("2e ints:\n$Ints")
    end
    for iter in 1:MaxIter
        D = dmat(U,nclosed)
        if verbose
            println("D=\n$D")
        end
        G = make2JmK(D,Ints)
        H = h+G
        E,U = eig(H,S)
        Eone = trace2(D,h)
        Etwo = trace2(D,H)
        Energy = Enuke + Eone + Etwo
        println("HF: $iter  $Energy : $Enuke    $Eone    $Etwo")
        if isapprox(Energy,Eold)
            break
        end
        Eold  = Energy
    end
    return Energy,E,U
end

function test_h2()
    @time Energy, E, U = rhf(h2)
    @assert isapprox(Energy,-1.1170996)
end

function test_lih()
    @time Energy, E, U = rhf(lih)
    @assert isapprox(Energy,-7.86073270525799)
end

function test_h2o()
    @time Energy,E,U = rhf(h2o)
    @assert isapprox(Energy,-74.9597609118851)
end

function test()
    test_utils()
    test_pgbf()
    test_cgbf()
    test_overlap()
    test_kinetic()
    test_a_terms()
    test_gamma()
    test_na()
    test_fgamma()
    test_one()
    test_na2()
    test_two_terms()
    test_coul1()
    test_vrr()
    test_hrr()
    test_geo_basis()
    test_h2()
    test_lih()
    test_h2o()
end

test()



