import dolfin as df
import math


def comp_axisymmetric_pure_neumann(mat_obj, mesh_obj, bc, omega, save_path):
    E = mat_obj.E 
    rho = mat_obj.rho 
    nu = mat_obj.nu 

    mesh = mesh_obj.create()
    Rext = mesh_obj.Rext
    Rint = mesh_obj.Rint

    # rename x[0], x[1] by x, y
    x, y = df.SpatialCoordinate(mesh)

    dim = mesh.topology().dim()


    cell_markers = df.MeshFunction("size_t", mesh, dim)
    facet_markers = df.MeshFunction("size_t", mesh, dim - 1)

    coord = mesh.coordinates()

    class OuterRadius(df.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0]**2 + x[1]**2 >= Rext**2 - 1e-3


    class InnerRadius(df.SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[0]**2 + x[1]**2 <= Rint**2 + 1e-3

    OuterRadius().mark(facet_markers, 1)
    InnerRadius().mark(facet_markers, 2)

    df.File(save_path + 'cell_markers.pvd') << cell_markers
    df.File(save_path + 'facet_markers.pvd') << facet_markers


    # Create mesh and define function space
    V = df.FunctionSpace(mesh, "CG", 1)
    # =========================== Define Finite Element ===========================
    degree = 1 #P1 table of periodic
    # CG: Continuous Galerkin
    V_ele = df.FiniteElement("CG", mesh.ufl_cell(), degree) # scalar, Ex: u
    R_ele = df.FiniteElement("R", mesh.ufl_cell(), 0) # scalar, Ex: u
    # Vector element # akternative : manually define vector function : better for 3D adaptation than above
    total_ele = df.MixedElement([V_ele, R_ele])
    W = df.FunctionSpace(mesh, total_ele)


    # Define trial function
    dunks = df.TrialFunction(W)
    (du, dc) = df.split(dunks)


    # Define variational problem
    tunks = df.TestFunction(W)
    (tu, tc) = df.split(tunks)


    # displacement in radial direction u(x,y)
    unks = df.Function(W)
    #(u, c) = df.split(unks)



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



    # Weak form
    dFu = r*sigma_r(du)*epsilon_r(tu)*df.dx(0)
    dFu = dFu + sigma_theta(du)*tu*df.dx(0)
    dFu = dFu + dc*tu*df.dx(0)
    dFu = dFu - rho*omega**2*r**2*tu*df.dx(0)
    dFc = du*tc*df.dx(1)



    dFu = r*sigma_r(du)*epsilon_r(tu)*df.dx
    dFu = dFu + sigma_theta(du)*tu*df.dx
    dFu = dFu + dc*tu*df.dx
    dFu = dFu - rho*omega**2*r**2*tu*df.dx
    dFc = du*tc*df.dx

    dF = dFu + dFc

    # residual 
    #F = df.action(dF, unks)


    a = df.lhs(dF)
    L = df.rhs(dF)


    df.solve(a == L, unks)

    (_u, _c) = unks.split(True)


    u_pro = df.project(_u, V)
    u_pro.rename('dis u [m]', 'dis u [m]')

    # displacement
    df.File(save_path + 'displacement.pvd') << u_pro


    # compute stresses
    sigma_r_pro = df.project(sigma_r(_u), V)
    sigma_r_pro.rename('sigma_r [Pa]', 'sigma_r [Pa]')
    df.File(save_path + 'sigma_r.pvd') << sigma_r_pro

    sigma_theta_pro = df.project(sigma_theta(_u), V) 
    sigma_theta_pro.rename('sigma_theta [Pa]', 'sigma_theta [Pa]')
    df.File(save_path + 'sigma_theta.pvd') << sigma_theta_pro

    # compute von Mises stress
    def von_mises_stress(sigma_r, sigma_theta):
        return df.sqrt(sigma_r**2 + sigma_theta**2 - sigma_r*sigma_theta) 

    von_stress_pro = df.project(von_mises_stress(sigma_r(_u), sigma_theta(_u)), V)
    von_stress_pro.rename('von Mises Stress [Pa]', 'von Mises Stress [Pa]')
    df.File(save_path + 'von_mises_stress.pvd') << von_stress_pro
















