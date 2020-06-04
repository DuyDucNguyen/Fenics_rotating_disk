import dolfin as df
import mshr 


def create_simple_mesh(Rint, Rout, resolution): 
    # Create circles as Circle(Center, Radius)
    outer_circle = mshr.Circle(df.Point(0,0), Rout)
    inner_circle = mshr.Circle(df.Point(0,0), Rint)
    domain = outer_circle - inner_circle
    mesh = mshr.generate_mesh(domain, resolution)
    return mesh 
