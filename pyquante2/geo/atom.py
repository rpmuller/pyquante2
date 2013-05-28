"""
Class to create an atom object

>>> h = atom(1,0,0,0)
>>> h
(1, 0.0, 0.0, 0.0)
>>> h.r
array([ 0.,  0.,  0.])

Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from numpy import array

class atom:
    def __init__(self,atno,x,y,z,atid=0,fx=0.0,fy=0.0,fz=0.0,vx=0.0,vy=0.0,vz=0.0,**kwargs):
        self.atno = atno
        self.r = array([x,y,z],'d')
        self.atid = atid
        self.f = array([fx,fy,fz],'d')
        self.vel = array([vx,vy,vz],'d')

        self.units = kwargs.get('units','bohr').lower()
        assert self.units[:4] in ['bohr','angs']
        if not self.units == 'bohr':
            self.r /= 0.52918

        return

    def atuple(self): return (self.atno,self.r[0],self.r[1],self.r[2])

    def __repr__(self): return repr(self.atuple())

    def __getitem__(self, i):
        return self.r[i]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
