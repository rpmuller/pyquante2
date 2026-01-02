#!/usr/bin/env python

import math
import numpy as np

from pyquante2.geo.elements import sym2no
from pyquante2.viewer.trackball_camera import TrackballCamera

try:
    from pyglet.gl import (
        glBegin, glEnd, glVertex3f, glColor3f, glPushMatrix, glPopMatrix,
        glTranslatef, glRotatef, glMatrixMode, glLoadIdentity, glViewport,
        glEnable, glDisable, glClear, glNewList, glEndList, glCallList,
        glLineWidth, glMaterialfv, glLightfv, glLightf,
        gluNewQuadric, gluQuadricDrawStyle, gluSphere, gluCylinder, gluPerspective, gluLookAt,
        GL_LINES, GL_LIGHTING, GL_LINE_SMOOTH, GL_DEPTH_TEST, GL_CULL_FACE, GL_PROJECTION,
        GL_MODELVIEW, GL_COMPILE, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT,
        GL_FRONT, GL_DIFFUSE, GL_SHININESS, GL_SPECULAR, GL_LIGHT0, GL_POSITION,
        GL_CONSTANT_ATTENUATION, GL_LINEAR_ATTENUATION, GLU_LINE, GLU_FILL,
        GLfloat, Config
    )
    from pyglet import window
except ImportError:
    pass


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

def draw_points(points):
    n,sb4 = points.shape
    pyglet.graphics.draw(n,pyglet.gl.GL_POINTS,
                         ('v3f',points.flatten()))
    return

def glf(x): return (GLfloat * len(x))(*x)

def norm1(x,maxx):
    """given x within [0,maxx], scale to a range [-1,1]."""
    return (2.0 * x - float(maxx)) / float(maxx)

class Viewer(object):
    def __init__(self,width=defaults['width'],height=defaults['height']):
        self.width = width
        self.height = height
        self.config = Config(double_buffer=True, depth_size=24)
        self.win = window.Window(visible=True,resizable=True,
                                 config=self.config, caption='PyQuante Viewer')

        # set callbacks. I think I could do all of these with
        # @window.event callbacks
        self.win.on_resize = self.on_resize
        self.win.on_draw = self.on_draw
        self.win.on_mouse_press = self.on_mouse_press
        self.win.on_mouse_drag = self.on_mouse_drag
        self.win.on_mouse_scroll = self.on_mouse_scroll

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
        from pyglet.window import key
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

    # This doesn't work at all
    def on_mouse_scroll(self,x,y,scroll_x, scroll_y):
        scalef=1.0
        dx = norm1(scalef*scroll_x,self.width)
        dy = norm1(scalef*scroll_y,self.height)
        #self.tb.mouse_zoom(dx,dy)
        return

class Points(object):
    def __init__(self,points): self.points = points
    def draw(self): draw_points(self.points)
        

class Sphere(object):
    def __init__(self,x,y,z,r,g,b,rad):
        self._pos = (x,y,z)
        self._color = (r,g,b)
        self._rad = rad

    def draw(self):
        x,y,z = self._pos
        r,g,b = self._color
        draw_sphere(x,y,z,r,g,b,self._rad)

    def __repr__(self):
        return "Sphere(%s,%s,%f)" % (self._pos,self._color,self._rad)

class Cylinder(object):
    def __init__(self,x1,y1,z1,x2,y2,z2,r,g,b,rad):
        self._start = (x1,y1,z1)
        self._end = (x2,y2,z2)
        self._color = (r,g,b)
        self._rad = rad

    def draw(self):
        x1,y1,z1 = self._start
        x2,y2,z2 = self._end
        r,g,b = self._color
        draw_cylinder(x1,y1,z1,x2,y2,z2,r,g,b,self._rad)

    def __repr__(self):
        return "Cylinder(%s,%s,%s,%f)" % (self._start,self._end,
                                          self._color,self._rad)

class Line(object):
    def __init__(self,x1,y1,z1,x2,y2,z2,r=255,g=255,b=255):
        self._start = (x1,y1,z1)
        self._end = (x2,y2,y2)
        self._color = (r,g,b)

    def draw(self):
        x1,y1,z1 = self._start
        x2,y2,z2 = self._end
        r,g,b = self._color
        draw_line(x1,y1,z1,x2,y2,z2,r,g,b)

class Shapes(object):
    def __init__(self,molecule=[],**kwargs):
        self.atoms = molecule
        self.shapelist = []
        for atom in self.atoms:
            self.add_atom(atom)
        self.find_bonds(**kwargs)
        for bond in self.bonds:
            self.add_bond(bond)
        return

    def add_sphere(self,x,y,z,r=0.9,g=0.9,b=0.9,rad=0.1):
        self.shapelist.append(Sphere(x,y,z,r,g,b,rad))

    def add_cylinder(self,x1,y1,z1,x2,y2,z2,r=0.5,g=0.5,b=0.5,rad=0.2):
        self.shapelist.append(Cylinder(x1,y1,z1,x2,y2,z2,r,g,b,rad))

    def add_points(self,points):
        #self.shapelist.append(Points(points)) # only visible from 1 side
        for i in range(points.shape[0]):
            x,y,z = points[i,:3]
            self.add_sphere(x,y,z,rad=0.02)
        return

    def add_points_weights(self,points):
        "Same as add_points, but color by 4th column value"
        from pyquante2.utils import colorscale
        weights = points[:,3]
        wmin = np.min(weights)
        wmax = np.max(weights)
        for i in range(points.shape[0]):
            x,y,z = points[i,:3]
            r,g,b = colorscale(weights[i],wmin,wmax)
            self.add_sphere(x,y,z,r,g,b,rad=0.02)
        return

    def add_atom(self,atom,**kwargs):
        x,y,z = atom.r
        r,g,b = atom.color()
        style = kwargs.get('style','BallStick')
        scaling = kwargs.get('scaling',0.5)
        rad = scaling*atom.radius()
        return self.add_sphere(x,y,z,r,g,b,rad)

    def add_bond(self,atoms,**kwargs):
        style = kwargs.get('style','BallStick')
        if style == 'Ball': return []
        x1,y1,z1 = atoms[0].r
        x2,y2,z2 = atoms[1].r
        r = kwargs.get('r',0.5)
        g = kwargs.get('g',0.5)
        b = kwargs.get('b',0.5)
        rad = kwargs.get('rad',0.2)
        return self.add_cylinder(x1,y1,z1,x2,y2,z2,r,g,b,rad)

    def find_bonds(self,scalef=0.6):
        from pyquante2.utils import upairs
        from pyquante2.constants import bohr2ang
        self.bonds = []
        for i,j in upairs(range(len(self.atoms))):
            ati,atj = self.atoms[i],self.atoms[j]
            r = ati.distance(atj)*bohr2ang
            r0 = ati.radius() + atj.radius()
            if r < scalef*r0:
                self.bonds.append((ati,atj))
        return

def test_prims():
    win = Viewer()
    spheres = [Sphere(-1,-1,0,1.,0.,0.,1.),
               Sphere(1,1,1.,0.,0.,1.,1.)]
    cyls = [Cylinder(-1,-1,0,1,1,1,0.5,0.5,0.5,0.2)]
    win.calllist(spheres+cyls)
    win.run()
    return

def test_mol():
    from pyquante2 import h2o
    shapes = Shapes(h2o)
    win = Viewer()
    win.calllist(shapes.shapelist)
    win.run()
    return


if __name__ == '__main__':
    #test_prims()
    test_mol()
