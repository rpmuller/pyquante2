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
    for ibf in range(len(bfs)):
        fxyz += orb[ibf]*bfs[ibf](x, y, z)
    fxyz = np.abs(fxyz)**2
    src = mlab.pipeline.scalar_field(x, y, z, fxyz)
    mlab.pipeline.iso_surface(src, contours=[fxyz.min()+0.02*fxyz.ptp(),], opacity=0.6)
    if doshow: mlab.show()
    return

def plot_orbs(molcule, orb, bfs):
    from mayavi import mlab
    #mlab.figure(1, bgcolor=(0, 0, 0), size=(750, 750))
    #mlab.clf()
    
    xarray = np.zeros((len(molecule), ))
    yarray = np.zeros((len(molecule), ))
    zarray = np.zeros((len(molecule), ))
    white = (1,1,1)
    gray = (0.5, 0.5, 0.5)
    red = (1, 0, 0)
    green = (0, 1, 0)
    blue = (0, 0, 1)
    color_dict = {1: (1, 1, 1), 7: blue, 6: gray, 8: red}
    scale_dict = {1: 1, 7: 1.5, 6: 1.5, 8: 1.5}

    for i in range(len(molecule)):
        (xarray[i], yarray[i], zarray[i]) = molecule[i].r
        at = mlab.points3d(xarray[i], yarray[i], zarray[i],
                       scale_factor=scale_dict[molecule[i].atno],
                       resolution=20,
                       color=color_dict[molecule[i].atno],
                       scale_mode='none')
    
    x, y, z = np.mgrid[min(xarray)-10.0:max(xarray)+10.0:100j,
                       min(yarray)-10.0:max(yarray)+10.0:100j, 
                       min(zarray)-10.0:max(zarray)+10.0:100j]
    fxyz = np.zeros((100, 100, 100))
    for ibf in range(len(bfs)):
        fxyz += orb[ibf]*bfs[ibf](x, y, z)
    fxyz = np.abs(fxyz)**2
    src = mlab.pipeline.scalar_field(x, y, z, fxyz)
    mlab.pipeline.iso_surface(src, contours=[fxyz.min()+0.02*fxyz.ptp(),], opacity=0.6)
    mlab.show()
