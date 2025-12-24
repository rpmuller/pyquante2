cdef extern from "cints.h":
    double overlap(double alpha1, int l1, int m1, int n1,
                   double xa, double ya, double za,
                   double alpha2, int l2, int m2, int n2,
                   double xb, double yb, double zb)
    double kinetic(double alpha1, int l1, int m1, int n1,
                   double xa, double ya, double za,
                   double alpha2, int l2, int m2, int n2,
                   double xb, double yb, double zb)
    double nuclear_attraction(double x1, double y1, double z1, double norm1,
                              int l1, int m1, int n1, double alpha1,
                              double x2, double y2, double z2, double norm2,
                              int l2, int m2, int n2, double alpha2,
                              double x3, double y3, double z3)
    double coulomb_repulsion(double xa, double ya, double za, double norma,
                             int la, int ma, int na, double alphaa,
                             double xb, double yb, double zb, double normb,
                             int lb, int mb, int nb, double alphab,
                             double xc, double yc, double zc, double normc,
                             int lc, int mc, int nc, double alphac,
                             double xd, double yd, double zd, double normd,
                             int ld, int md, int nd, double alphad)


