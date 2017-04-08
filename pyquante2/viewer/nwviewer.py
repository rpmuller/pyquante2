from pyquante2.viewer.viewer import Viewer

class NWViewer(Viewer):
    def __init__(self):
        Viewer.__init__(self)
        return

def run():
    win = NWViewer()
    # Consider parsing args here
    win.run()
    return

if __name__ == '__main__': run()


