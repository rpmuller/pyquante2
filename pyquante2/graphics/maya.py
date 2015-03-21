import numpy as np

def view_dft_density(grid,D,bbox,npts=50,doshow=True):
    from mayavi import mlab
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


# Thanks to Thomas Markovich for this code:
def view_mol(mol,doshow=True):
    from pyquante2.element import color,radius
    from mayavi import mlab
    for at in mol:
        rgb = tuple(c/255. for c in color[at.Z])
        mlab.points3d(at.r[0],at.r[1],at.r[2],
                      scale_factor=radius[at.Z],color=rgb,
                      resolution=20,scale_mode='none')
    # Draw bonds?
    # Draw in cylinder mode?
    if doshow: mlab.show()
    return

def view_orb(mol,orb,bfs,npts=50,doshow=True):
    from mayavi import mlab
    xmin,xmax,ymin,ymax,zmin,zmax = mol.bbox()
    x, y, z = np.mgrid[xmin:xmax:(npts*1j),ymin:ymax:(npts*1j),zmin:zmax:(npts*1j)]

    fxyz = np.zeros((npts, npts, npts))
    for i,bf in enumerate(bfs):
        fxyz += orb[i]*bf(x, y, z)
    fxyz = np.abs(fxyz)**2
    src = mlab.pipeline.scalar_field(x, y, z, fxyz)
    mlab.pipeline.iso_surface(src, contours=[fxyz.min()+0.02*fxyz.ptp(),], opacity=0.6)
    if doshow: mlab.show()
    return

