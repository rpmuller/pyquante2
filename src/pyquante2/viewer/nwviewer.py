#!/usr/bin/env python
"""Viewer/editor for nwchem files.
"""

from pyquante2.viewer.viewer import Viewer

# For 2-3 compatibility:
try:
   input = raw_input
except NameError:
   pass

class NWViewer(Viewer):
    def __init__(self):
        Viewer.__init__(self)
        self.win.on_key_press = self.on_key_press
        return

    def on_key_press(self, key, modifiers):
        if key == window.key.I:
            self.input_geo()
        elif key == window.key.LEFT:
            self.prev_geo()
        elif key == window.key.RIGHT:
            self.next_geo()
        elif key == window.key.Q:
            self.quit()
        return

    def quit(self): pyglet.app.exit()

    def input_geo(self):
        fname = input("filename: ")
        print(fname)
        return

def run():
    win = NWViewer()
    # Consider parsing args here
    win.run()
    return

if __name__ == '__main__': run()


