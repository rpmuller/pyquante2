from pyglet import gl
# Defaults for kwarg-passed options
defaults = {
    'width'       : 600,  # Window width
    'height'      : 600,  # Window height
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
    'lightpos'    : (10,4,10,1), # Default position for lighting
    'lightcolor'  : (1,1,1,1),   # Default color for lighting
    }

def glf(x): return (gl.GLfloat * len(x))(*x)

def norm1(x,maxx):
    """given x within [0,maxx], scale to a range [-1,1]."""
    return (2.0 * x - float(maxx)) / float(maxx)

