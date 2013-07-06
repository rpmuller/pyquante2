#!/usr/bin/env python

from pyglet.gl import *
from pyglet import window
from trackball_camera import TrackballCamera
from Element import Element,sym2no
import math
import numpy

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

def parseline(line,format):
    """\
    Given a line (a string actually) and a short string telling
    how to format it, return a list of python objects that result.

    The format string maps words (as split by line.split()) into
    python code:
    x   ->    Nothing; skip this word
    s   ->    Return this word as a string
    i   ->    Return this word as an int
    d   ->    Return this word as an int
    f   ->    Return this word as a float

    Basic parsing of strings:
    >>> parseline('Hello, World','ss')
    ['Hello,', 'World']

    You can use 'x' to skip a record; you also don't have to parse
    every record:
    >>> parseline('1 2 3 4','xdd')
    [2, 3]

    >>> parseline('C1   0.0  0.0 0.0','sfff')
    ['C1', 0.0, 0.0, 0.0]

    Should this return an empty list?
    >>> parseline('This line wont be parsed','xx')
    """
    xlat = {'x':None,'s':str,'f':float,'d':int,'i':int}
    result = []
    words = line.split()
    for i in range(len(format)):
        f = format[i]
        trans = xlat.get(f,None)
        if trans: result.append(trans(words[i]))
    if len(result) == 0: return None
    if len(result) == 1: return result[0]
    return result

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

class TBWindow:
    def __init__(self,width=defaults['width'],height=defaults['height']):
        self.width = width
        self.height = height
        self.config = Config(double_buffer=True, depth_size=24)
        self.win = window.Window(visible=True,resizable=True,
                                 config=self.config, caption='Sam')

        # set callbacks
        self.win.on_resize = self.on_resize
        self.win.on_draw = self.on_draw
        self.win.on_mouse_press = self.on_mouse_press
        self.win.on_mouse_drag = self.on_mouse_drag

        self.win.set_size(self.width,self.height)

        self.init_gl()
        
        self.tb = TrackballCamera(20.0)
        self.clnum = 1
        return

    def init_gl(self,**kwargs):
        swire_draw = kwargs.get('swire_draw',defaults['swire_draw'])
        mat_draw = kwargs.get('mat_draw',defaults['mat_draw'])
        glEnable(GL_DEPTH_TEST)
        glDisable(GL_CULL_FACE)
        if swire_draw:
            glEnable(GL_LINE_SMOOTH)
        elif mat_draw:
            glEnable(GL_LIGHTING)
            lightZeroPosition = glf(defaults['lightpos'])
            lightZeroColor = glf(defaults['lightcolor'])
            glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
            glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
            glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
            glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
            glEnable(GL_LIGHT0)
        return

    def run(self): pyglet.app.run()

    def calllist(self, shapes):
        glNewList(self.clnum,GL_COMPILE)
        for shape in shapes: shape.draw()
        glEndList()
        return

    def on_resize(self, width, height):
        self.width = width
        self.height = height
        glViewport(0,0,self.width,self.height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective( 
            40.0,                            # Field Of View
            float(self.width)/float(self.height),  # aspect ratio
            1.0,                             # z near
            100.0)                           # z far
        self.tb.update_modelview() # init modview matrix for trackball
        return

    def on_draw(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glCallList(self.clnum)
        return

    def on_mouse_press(self, x, y, button, modifiers):
        if button == window.mouse.LEFT:
            self.tb.mouse_roll(
                norm1(x, self.width),
                norm1(y,self.height),
                False)
        elif button == window.mouse.RIGHT:
            self.tb.mouse_zoom(
                norm1(x, self.width),
                norm1(y,self.height),
                False)
        return

    def on_mouse_drag(self, x, y, dx, dy, buttons, modifiers):
        if buttons & window.mouse.LEFT:
            self.tb.mouse_roll(
                norm1(x,self.width),
                norm1(y,self.height))
        elif buttons & window.mouse.RIGHT:
            self.tb.mouse_zoom(
                norm1(x,self.width),
                norm1(y,self.height))
        return

class Position:
    def __init__(self,x=0,y=0,z=0):
        self._pos = numpy.array((x,y,z),'d')
        return

    def get(self): return self._pos

    def translate(self,d):
        if type(d) == type((0,)):
            self._pos += numpy.array(d,'d')
        self._pos += d
        return

    def dist(self,other):
        d = self._pos - other.get()
        return math.sqrt(numpy.dot(d,d))

    def __repr__(self): return "(%.4f,%.4f,%.4f)" % tuple(self._pos)

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
    
class Material:
    def __init__(self):
        self._atoms = []
        self._bonds = []
        self._uc = None

    def nat(self): return len(self._atoms)

    def atoms(self,i=None):
        if i is not None: return self._atoms.__getitem__(i)
        return self._atoms
    def addatom(self,atom): self._atoms.append(atom)
    def bonds(self): return self._bonds
    def uc(self): return self._uc
    def translate(self,delta):
        for atom in self._atoms:
            atom.translate(delta)
        return

    def center(self):
        com = self.com()
        self.translate(-com)
        return

    def shapes(self,**kwargs):
        s = []
        for atom in self.atoms():
            s.extend(atom.shapes(**kwargs))
        for bond in self.bonds():
            s.extend(bond.shapes(**kwargs))
        if self._uc:
            s.extend(self._uc.shapes(**kwargs))
        return s

    def find_bonds(self,scalef=0.6):
        self._bonds = []
        for i,j in pairs(self.nat()):
            r = self.atoms(i).dist(self.atoms(j))
            r0 = self.atoms(i).radius() + self.atoms(j).radius()
            if r < scalef*r0:
                self._bonds.append(Bond(self.atoms(i),self.atoms(j)))
        return

    def com(self):
        # This is ugly: I should be working with positions, rather
        # than the underlying ndarrays
        p = Position()
        p = p._pos
        totm = 0
        for atom in self._atoms:
            totm += atom.mass()
            p += atom.mass()*atom.pos()
        return p/totm

def pairs(n,diag=False):
    d = 0
    if diag: d = 1
    for i in xrange(n):
        for j in xrange(i+d):
            yield i,j
    return

class Atom:
    def __init__(self,x,y,z,atno):
        self._type = Element(atno)
        self._pos = Position(x,y,z)
        self.ball_stick_scale = 0.3
        self.ball_scale = 0.9

    def dist(self,other): return self._pos.dist(other._pos)
    def radius(self): return self._type.radius()
    def color(self): return self._type.color()
    def pos(self): return self._pos.get()
    def mass(self): return self._type.mass()

    def scale(self,style):
        if style == 'Ball':
            return self.ball_scale
        return self.ball_stick_scale

    def translate(self,delta): self._pos.translate(delta)

    def shapes(self,**kwargs):
        x,y,z = self.pos()
        r,g,b = self.color()
        style = kwargs.get('style','BallStick')
        scaling = self.scale(style)
        rad = scaling*self.radius()
        return [Sphere(x,y,z,r,g,b,rad)]

class Bond:
    def __init__(self,start,end):
        self._start = start
        self._end = end
        return

    def shapes(self,**kwargs):
        style = kwargs.get('style','BallStick')
        if style == 'Ball': return []
        x1,y1,z1 = self._start.pos()
        x2,y2,z2 = self._end.pos()
        r = kwargs.get('r',0.5)
        g = kwargs.get('g',0.5)
        b = kwargs.get('b',0.5)
        rad = kwargs.get('rad',0.2)
        return [Cylinder(x1,y1,z1,x2,y2,z2,r,g,b,rad)]

def test_prims():
    win = TBWindow()
    spheres = [Sphere(-1,-1,0,1.,0.,0.,1.),
               Sphere(1,1,1.,0.,0.,1.,1.)]
    cyls = [Cylinder(-1,-1,0,1,1,1,0.5,0.5,0.5,0.2)]
    win.calllist(spheres+cyls)
    win.run()
    return

def read_xyz(fname):
    """
    Read an xyz format file into the sam structure.
    Currently returns only the final geometry.
    """
    f = open(fname)
    while True:
        line = f.readline()
        if not line: break
        nat = parseline(line,'d')
        mat = Material()
        comments = f.readline()
        for i in range(nat):
            line = f.readline()
            sym,x,y,z = parseline(line,'sfff')
            atno = sym2no[sym]
            mat.addatom(Atom(x,y,z,atno))
    mat.find_bonds()
    mat.center()
    return mat

def test_read():
    aspirin = "../Geos/aspirin.xyz"
    benzene = "../Geos/benzene.xyz"
    caffeine = "../Geos/caffeine.xyz"
    win = TBWindow()
    mat = read_xyz(aspirin)
    shapes = mat.shapes()
    win.calllist(shapes)
    win.run()
    return

def main():
    #test_prims()
    test_read()
    return

if __name__=="__main__": main()
