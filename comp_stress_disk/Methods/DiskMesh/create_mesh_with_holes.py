import dolfin as df
import mshr 


def create_mesh_with_holes(Rint, Rout, Rhole, resolution):
    # Create circles as Circle(Center, Radius)
    outer_circle = mshr.Circle(df.Point(0,0), Rout)
    inner_circle = mshr.Circle(df.Point(0,0), Rint)
    hole = mshr.Circle(df.Point(0, 0.5*(Rint+Rout)), Rhole)
    domain = outer_circle - inner_circle - hole
    mesh = mshr.generate_mesh(domain, resolution)
    return mesh 
