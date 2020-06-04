import dolfin as df
import math


#save_path = "comp_no_t/"

def define_markers(mesh, Rext, Rint):
    '''
    define markers
    
    :param mesh: mesh
    :param Rext: exterior radius
    :param Rint: inner radius
    '''
    dim = mesh.topology().dim()

    cell_markers = df.MeshFunction("size_t", mesh, dim)
    facet_markers = df.MeshFunction("size_t", mesh, dim - 1)

    #coord = mesh.coordinates()

    class OuterRadius(df.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0]**2 + x[1]**2 >= Rext**2 - 1e-3


    class InnerRadius(df.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0]**2 + x[1]**2 <= Rint**2 + 1e-3

    OuterRadius().mark(facet_markers, 1)
    InnerRadius().mark(facet_markers, 2)

    #df.File(save_path + 'facet_markers.pvd') << facet_markers

    return cell_markers, facet_markers

















