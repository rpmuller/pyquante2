import pyglet
from pyglet.gl import *

# Defaults for kwarg-passed options
defaults = {
    'mat_draw'    : True, # Use materials/lighting to draw
    'wire_draw'   : False,# Use open wires to draw
    'swire_draw'  : False,# Use smoothed wires to draw
    'fwire_draw'  : False,# Use filled wires to draw
    'specular'    : True, # Use specular highlights
    'sslices'     : 32,   # Number of slices in a sphere
    'sstacks'     : 32,   # Number of stacks in a sphere
    'cslices'     : 12,   # Number of slices in a cylinder
    'cstacks'     : 12,   # Number of stacks in a cylinder
    'lweight'     : 1,    # Line weight
    'width'       : 600,  # Window width
    'height'      : 600,  # Window height
    'lightpos'    : (10,4,10,1), # Default position for lighting
    'lightcolor'  : (1,1,1,1),   # Default color for lighting
    }


def draw_sphere(x,y,z,red,green,blue,rad,**kwargs):
    mat_draw = kwargs.get('mat_draw',defaults['mat_draw'])
    fwire_draw = kwargs.get('fwire_draw',defaults['fwire_draw'])
    sslices = kwargs.get('sslices',defaults['sslices'])
    sstacks = kwargs.get('sstacks',defaults['sstacks'])
    specular = kwargs.get('specular',defaults['specular'])
    glColor3f(red, green, blue)
    q = gluNewQuadric()
    if mat_draw:
        glMaterialfv(GL_FRONT, GL_DIFFUSE, glf((red,green,blue,1.0)))
        if specular:
            glMaterialfv(GL_FRONT, GL_SHININESS, glf([25.0]))
            glMaterialfv(GL_FRONT, GL_SPECULAR, glf([1.0,1.0,1.0,1.0])) 
    elif fwire_draw:
        gluQuadricDrawStyle(q,GLU_FILL)
    else:
        gluQuadricDrawStyle(q,GLU_LINE)
    glPushMatrix()
    glTranslatef(x, y, z)
    gluSphere(q,rad,sslices,sstacks)
    glPopMatrix()
    return

def draw_grid():
    # Either remove hardwires or delete
    glBegin(GL_LINES)
    glColor3f(1.0, 1.0, 1.0)
    for i in range(10):
        glVertex3f(i*10.0,-100., 0.)
        glVertex3f(i*10.0, 100., 0.)
        
        glVertex3f(-i*10.0,-100., 0.)
        glVertex3f(-i*10.0, 100., 0.)

        glVertex3f(-100., i*10.0, 0.)
        glVertex3f( 100., i*10.0, 0.)

        glVertex3f(-100.,-i*10.0, 0.)
        glVertex3f( 100.,-i*10.0, 0.)
    glEnd()
    return

def draw_cylinder(x1,y1,z1,x2,y2,z2,red,green,blue,rad,**kwargs):
    mat_draw = kwargs.get('mat_draw',defaults['mat_draw'])
    fwire_draw = kwargs.get('fwire_draw',defaults['fwire_draw'])
    cslices = kwargs.get('cslices',defaults['cslices'])
    cstacks = kwargs.get('cstacks',defaults['cstacks'])
    specular = kwargs.get('specular',defaults['specular'])

    rad2deg=180.0/math.pi
    dx,dy,dz = x2-x1,y2-y1,z2-z1
    length = math.sqrt(dx*dx+dy*dy+dz*dz)
    dxn,dyn,dzn = dx/length,dy/length,dz/length
    rx,ry,rz = -dyn,dxn,0
    theta = math.acos(dzn)*rad2deg

    glPushMatrix()
    glColor3f(red, green, blue)
    q = gluNewQuadric()
    if mat_draw:
        glMaterialfv(GL_FRONT, GL_DIFFUSE, glf((red,green,blue,1.0)))
        if specular:
            glMaterialfv(GL_FRONT, GL_SHININESS, glf([25.0]))
            glMaterialfv(GL_FRONT, GL_SPECULAR, glf([1.0,1.0,1.0,1.0])) 
    elif fwire_draw:
        gluQuadricDrawStyle(q,GLU_FILL)
    else:
        gluQuadricDrawStyle(q,GLU_LINE)
    glTranslatef(x1,y1,z1)

    glRotatef(theta,rx,ry,rz)
    gluCylinder(q,rad,rad,length,cslices,cstacks)
    glPopMatrix()
    return

def draw_line(x1,y1,z1,x2,y2,z2,red,green,blue,**kwargs):
    lweight = kwargs.get('lweight',defaults['lweight'])
    glDisable(GL_LIGHTING)
    glEnable(GL_LINE_SMOOTH)
    glLineWidth(lweight)
    glColor3f(red,green,blue)
    glBegin(GL_LINES)
    glVertex3f(x1,y1,z1)
    glVertex3f(x2,y2,z2)
    glEnd()
    glDisable(GL_LINE_SMOOTH)
    glEnable(GL_LIGHTING)
    return

def glf(x): return (GLfloat * len(x))(*x)

def norm1(x,maxx):
    """given x within [0,maxx], scale to a range [-1,1]."""
    return (2.0 * x - float(maxx)) / float(maxx)


class Sphere:
    def __init__(self,x,y,z,r,g,b,rad):
        self._pos = Position(x,y,z)
        self._color = Color(r,g,b)
        self._rad = rad

    def draw(self):
        x,y,z = self._pos.get()
        r,g,b = self._color.get()
        draw_sphere(x,y,z,r,g,b,self._rad)

    def __repr__(self):
        return "Sphere(%s,%s,%f)" % (self._pos,self._color,self._rad)

class Color:
    def __init__(self,r,g,b):
        self._color = (r,g,b)
        return

    def get(self): return self._color

    def __repr__(self): return "(%.2f,%.2f,%.2f)" % self._color

class Cylinder:
    def __init__(self,x1,y1,z1,x2,y2,z2,r,g,b,rad):
        self._start = Position(x1,y1,z1)
        self._end = Position(x2,y2,z2)
        self._color = Color(r,g,b)
        self._rad = rad

    def draw(self):
        x1,y1,z1 = self._start.get()
        x2,y2,z2 = self._end.get()
        r,g,b = self._color.get()
        draw_cylinder(x1,y1,z1,x2,y2,z2,r,g,b,self._rad)

    def __repr__(self):
        return "Cylinder(%s,%s,%s,%f)" % (self._start,self._end,
                                          self._color,self._rad)

class Line:
    def __init__(self,x1,y1,z1,x2,y2,z2,r=255,g=255,b=255):
        self._start = Position(x1,y1,z1)
        self._end = Position(x2,y2,y2)
        self._color = Color(r,g,b)

    def draw(self):
        x1,y1,z1 = self._start.get()
        x2,y2,z2 = self._end.get()
        r,g,b = self._color.get()
        draw_line(x1,y1,z1,x2,y2,z2,r,g,b)

class UC:
    def __init__(self,A,B,C,origin=(0,0,0)):
        self._origin = Position(*origin)
        self._A = Position(*A)
        self._B = Position(*B)
        self._C = Position(*C)
        return

    def shapes(self):
        xo,yo,zo = self._origin.get()
        xa,ya,za = self._A.get()
        xb,yb,zb = self._B.get()
        xc,yc,zc = self._C.get()
        s = [Line(xo,yo,zo,xo+xa,yo+ya,zo+za),
             Line(xo,yo,zo,xo+xb,yo+yb,zo+zb),
             Line(xo,yo,zo,xo+xc,yo+yc,zo+zc),
             Line(xo+xa+xb+xc,zo+ya+yb+yc,zo+za+zb+zc,
                  xo+xa,yo+ya,zo+za),
             Line(xo+xa+xb+xc,zo+ya+yb+yc,zo+za+zb+zc,
                  xo+xb,yo+yb,zo+zb),
             Line(xo+xa+xb+xc,zo+ya+yb+yc,zo+za+zb+zc,
                  xo+xc,yo+yc,zo+zc)]
        return s
    
class Shapes:
    def __init__(self,molecule,**kwargs):
        self.atoms = molecule
        self.find_bonds(**kwargs)
        return

    def shapes(self,**kwargs):
        s = []
        for atom in self.atoms():
            s.extend(self.atom_shapes(atom,**kwargs))
        for bond in self.bonds():
            s.extend(self.bond_shapes(bond,**kwargs))
        return s

    def atom_shapes(self,atom,**kwargs):
        x,y,z = self.pos()
        r,g,b = self.color()
        style = kwargs.get('style','BallStick')
        scaling = self.scale(style)
        rad = scaling*self.radius()
        return [Sphere(x,y,z,r,g,b,rad)]

    def bond_shapes(self,bond,**kwargs):
        style = kwargs.get('style','BallStick')
        if style == 'Ball': return []
        x1,y1,z1 = self._start.pos()
        x2,y2,z2 = self._end.pos()
        r = kwargs.get('r',0.5)
        g = kwargs.get('g',0.5)
        b = kwargs.get('b',0.5)
        rad = kwargs.get('rad',0.2)
        return [Cylinder(x1,y1,z1,x2,y2,z2,r,g,b,rad)]

    def find_bonds(self,scalef=0.6):
        from pyquante2.utils import upairs
        self.bonds = []
        for i,j in upairs(len(self.atoms)):
            ati,atj = self.atoms[i],self.atoms[j]
            r = ati.distance(atj)
            r0 = ati.radius() + atj.radius()
            if r < scalef*r0:
                self.bonds.append(Bond(self.atoms(i),self.atoms(j)))
        return


def test_prims():
    win = TBWindow()
    spheres = [Sphere(-1,-1,0,1.,0.,0.,1.),
               Sphere(1,1,1.,0.,0.,1.,1.)]
    cyls = [Cylinder(-1,-1,0,1,1,1,0.5,0.5,0.5,0.2)]
    win.calllist(spheres+cyls)
    win.run()
    return


if __name__ == '__main__':
    viewer()

