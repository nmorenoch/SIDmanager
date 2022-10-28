"""

@author: Nicolas Moreno Chaparro
         nicolas.morenochaparro@kaust.edu.sa

It generates 3d contour plots usign mayavi funcitonalities.
input: Volumetric data (3d array), plane show a slice, and position of the slice

"""


from mayavi import mlab



def plot3D(F,  planeO, sI = 0):
    mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
    mlab.clf()

        
    #mlab.contour3d(O, contours=3)
    #mlab.show()
        
    src = mlab.pipeline.scalar_field(F)
    min = F.min()
    max = F.max()
    
    mlab.pipeline.iso_surface(src, contours=[F.min()+0.8*(max-min), ], opacity=0.1)
    #mlab.pipeline.iso_surface(src, contours=[F.min()+0.1*F.ptp(), ], opacity=0.1)
    mlab.pipeline.iso_surface(src, contours=[F.max(), ],)
    mlab.pipeline.image_plane_widget(src,
                           plane_orientation=planeO,
                            slice_index=sI,
                        )
        
    # vol = mlab.pipeline.volume(source, vmin=min+0.65*(max-min),
    #                            vmax=min+0.9*(max-min))
    #mlab.view(132, 54, 45, [21, 20, 21.5])
    # mlab.show()
    
    #pts = mlab.points3d(self.glob[:,0], self.glob[:,1], self.glob[:,2], 1.5*self.glob[:,3].max() - self.glob[:,3],
    #         scale_factor=0.015, resolution=10)
                                #vol = mlab.pipeline.volume(mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts)))
    mlab.show()
        
    #m = vs.VolumeSlicer(data=O)
    #m.configure_traits()
