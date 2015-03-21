import numpy as np
from mayavi import mlab

def view_dft_density(grid,D,bbox,npts=50,doshow=True):
    rho = grid.getdens_interpolated(D,bbox,npts)
    # mlab.contour3d(rho,contours=8,opacity=0.5)
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
                                     plane_orientation='x_axes',
                                     slice_index=25)
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
                                     plane_orientation='y_axes',
                                     slice_index=25)
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
                                     plane_orientation='z_axes',
                                     slice_index=25)
    if doshow: mlab.show()
    return
