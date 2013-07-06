
def test_mesh():
    from pyquante2 import grid,h2o
    from pyquante2.viewer.viewer import Shapes,Viewer
    h2o_mesh = grid(h2o)
    shapes = Shapes(h2o)
    shapes.add_points(h2o_mesh.points[:,:3])
    win = Viewer()
    win.calllist(shapes.shapelist)
    win.run()
    return

if __name__ == '__main__': test_mesh()
    
