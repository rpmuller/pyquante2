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
def view_mol(mol,doshow=True,scale=1.0):
    from mayavi import mlab
    for at in mol:
        mlab.points3d(at.r[0],at.r[1],at.r[2],
                      scale_factor=scale*at.radius(),color=at.color(),
                      resolution=20,scale_mode='none')
    # Draw bonds?
    # Draw in cylinder mode?
    if doshow: mlab.show()
    return

def view_bonds(mol,doshow=True):
    from mayavi import mlab
    for i,j in mol.bonds():
        mlab.plot3d([mol.atoms[i].r[0], mol.atoms[j].r[0]],
                    [mol.atoms[i].r[1], mol.atoms[j].r[1]],
                    [mol.atoms[i].r[2], mol.atoms[j].r[2]],
                    tube_radius=0.2,
                    tube_sides=20)
    if doshow: mlab.show()
    return

def view_orb(mol,orb,bfs,npts=50,posval=0.05,doshow=True,planes=[]):
    from mayavi import mlab
    xmin,xmax,ymin,ymax,zmin,zmax = mol.bbox()
    x, y, z = np.mgrid[xmin:xmax:(npts*1j),ymin:ymax:(npts*1j),zmin:zmax:(npts*1j)]

    fxyz = np.zeros((npts, npts, npts))
    for c,bf in zip(orb,bfs):
        fxyz += c*bf(x, y, z)

    src = mlab.pipeline.scalar_field(x, y, z, fxyz)
    mlab.pipeline.iso_surface(src, contours=[-posval,posval], opacity=0.6)
    for d,ind in planes:
        if not ind: ind = npts//2
        mlab.pipeline.image_plane_widget(src,
                                         plane_orientation='%s_axes' % d,
                                         slice_index=ind)
        
    if doshow: mlab.show()
    return

