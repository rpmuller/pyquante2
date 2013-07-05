from pyglet import *
from trackball_camera import TrackballCamera
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
