import dolfin as df
import math
import matplotlib.pyplot as plt

from .define_markers import define_markers

# Implement the model in Mohammad2017
# cartesian coordinate (x, y) 

def comp_axisymmetric_dirichlet(mat_obj, mesh_obj, bc, omega, save_path):
    E = mat_obj.E
    rho = mat_obj.rho
    nu = mat_obj.nu

    mesh = mesh_obj.create()
    Rext = mesh_obj.Rext
    Rint = mesh_obj.Rint

    cell_markers, facet_markers = define_markers(mesh, Rext, Rint)


    # rename x[0], x[1] by x, y
    x, y = df.SpatialCoordinate(mesh)

    dim = mesh.topology().dim()

    coord = mesh.coordinates()


    # Create mesh and define function space
    V = df.FunctionSpace(mesh, "CG", 1)

    # Define boundary condition (homogeneous BC)
    u0 = df.Constant(0.0)

    if bc == 'CC':
        bc = [df.DirichletBC(V, u0, facet_markers, i) for i in [1, 2]]
    elif bc == 'CF':
        bc = df.DirichletBC(V, u0, facet_markers, 2)


    # Define variational problem
    du = df.TrialFunction(V)
    tu = df.TestFunction(V)
    # displacement in radial direction u(x,y)
    u = df.Function(V, name='displacement')

    class THETA(df.UserExpression):
        def eval(self, values, x):
            values[0] = math.atan2(x[1],x[0])
        def value_shape(self):
            #return (1,) # vector
            return () # scalar

    theta = THETA(degree = 1)
    #theta_int = df.interpolate(theta, df.FunctionSpace(mesh, "DG", 0))
    #df.File(save_path + 'theta.pvd') << theta_int

    class RADIUS(df.UserExpression):
        def eval(self, values, x):
            values[0] = df.sqrt(x[0]*x[0]+x[1]*x[1])
        def value_shape(self):
            return () # scalar
    r = RADIUS(degree = 1)



    # strain radial
    def epsilon_r(du):
        return 1.0/r*(x*df.Dx(du, 0) + y*df.Dx(du, 1))

    # strain circumferential
    def epsilon_theta(du): 
        return du/df.sqrt(x**2 + y**2)

    # radial stress # train-stress relation
    def sigma_r(du):
        return E/(1.0-nu**2)*(epsilon_r(du) + nu*epsilon_theta(du))

    # circumferential stress
    def sigma_theta(du): 
        return E/(1.0-nu**2)*(nu*epsilon_r(du) + epsilon_theta(du))

    #define centrifugal force vector form
    Fc = rho*omega**2*r**2


    # Weak form
    dres = r*sigma_r(du)*epsilon_r(tu)*df.dx
    dres = dres + sigma_theta(du)*tu*df.dx
    dres = dres - Fc*tu*df.dx

    # residual 
    res = df.action(dres, u)

    # solve 
    df.solve(res == 0, u, bc)

    # displacement
    df.File(save_path + 'displacement.pvd') << u


    # compute stresses
    sigma_r_pro = df.project(sigma_r(u), V)
    sigma_r_pro.rename('sigma_r [Pa]', 'sigma_r [Pa]')
    df.File(save_path + 'sigma_r.pvd') << sigma_r_pro

    sigma_theta_pro = df.project(sigma_theta(u), V) 
    sigma_theta_pro.rename('sigma_theta [Pa]', 'sigma_theta [Pa]')
    df.File(save_path + 'sigma_theta.pvd') << sigma_theta_pro

    # compute von Mises stress
    def von_mises_stress(sigma_r, sigma_theta):
        return df.sqrt(sigma_r**2 + sigma_theta**2 - sigma_r*sigma_theta) 

    von_stress_pro = df.project(von_mises_stress(sigma_r(u), sigma_theta(u)), V)
    von_stress_pro.rename('von Mises Stress [Pa]', 'von Mises Stress [Pa]')
    df.File(save_path + 'von_mises_stress.pvd') << von_stress_pro


    # save results to h5
    rfile = df.HDF5File(mesh.mpi_comm(), save_path + 'results.h5', "w")
    rfile.write(u, "u")
    rfile.write(sigma_r_pro, "sigma_r")
    rfile.write(sigma_theta_pro, "sigma_theta")
    rfile.write(von_stress_pro, "von_mises_stress")
    rfile.close()


    '''
    # Load solution
    U = df.Function(V)
    Sig_r = df.Function(V)
    Sig_theta = df.Function(V)
    Vm_stress = df.Function(V)
    input_file = df.HDF5File(mesh.mpi_comm(), save_path + 'results.h5', "r")
    input_file.read(U, "u")
    input_file.read(Sig_r, "sigma_r")
    input_file.read(Sig_theta, "sigma_theta")
    input_file.read(Vm_stress, "von_mises_stress")
    input_file.close()
    '''

    return 

















