import dolfin as df
import math
import matplotlib.pyplot as plt
from .define_markers import define_markers
from ufl import *

# Vullo - Rotors: Stess Analysis and Design
# Full model of u, v and thermal load T
# cartesian coordinate (x, y) 



def comp_full_model_heat_dirichlet(mat_obj, mesh_obj, bc, omega, save_path):
    # ======
    # Parameters
    # ======
    E = mat_obj.E
    rho = mat_obj.rho
    nu = mat_obj.nu

    mesh = mesh_obj.create()
    Rext = mesh_obj.Rext
    Rint = mesh_obj.Rint
    G = 1 # FAKE


    # ======
    # Thermal load
    # ======
    alpha = 1
    T = 1


    # ======
    # Thickness profile
    # ======
    h = 1 


    omega_velo = 1 #Fake


    # ======
    # markers
    # ======
    cell_markers, facet_markers = define_markers(mesh, Rext, Rint)

    
    # rename x[0], x[1] by x, y
    x, y = df.SpatialCoordinate(mesh)

    dim = mesh.topology().dim()

    coord = mesh.coordinates()

    # ======
    # Create function space
    # ======
    V = df.FunctionSpace(mesh, "CG", 1)
    degree = 1 
    fi_ele = FiniteElement("CG", mesh.ufl_cell(), degree) # CG: Continuous Galerkin
    vec_ele = VectorElement("CG", mesh.ufl_cell(), degree) 
    total_ele = MixedElement([fi_ele, fi_ele])
    W = df.FunctionSpace(mesh, total_ele)

  

    # ======
    # Define boundary condition
    # ======
    if bc == 'CC':
        u_Dbc = [df.DirichletBC(W.sub(0), df.Constant(0.0), facet_markers, bc) for bc in [1, 2]]
        v_Dbc = [df.DirichletBC(W.sub(1), df.Constant(0.0), facet_markers, bc) for bc in [1, 2]]
    elif bc == 'CF':
        u_Dbc = [df.DirichletBC(W.sub(0), df.Constant(0.0), facet_markers, bc) for bc in [2]]
        v_Dbc = [df.DirichletBC(W.sub(1), df.Constant(0.0), facet_markers, bc) for bc in [2]]
    Dbc = u_Dbc + v_Dbc



    # ======
    # Define functions
    # ======
    dunks = df.TrialFunction(W)
    tunks = df.TestFunction(W)
    unks = df.Function(W, name='displacement')

    # u(x,y): displacement in radial direction 
    # v(x,y): displacement in tangential direction 
    (du, dv) = df.split(dunks)
    (tu, tv) = df.split(tunks)
    (u, v) = df.split(unks)

    

    # ======
    # Define variable
    # ======
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


    # ======
    # Define week form
    # ======
    def d_dr(du):
        return 1.0/r*(x*df.Dx(du, 0) + y*df.Dx(du, 1))

    def d_dtheta(du):
        return -y*df.Dx(du, 0) + x*df.Dx(du, 1)

    # strain radial
    def epsilon_r(du):
        return d_dr(du)

    # strain circumferential
    def epsilon_theta(du, dv): 
        return du/r + 1.0/r*d_dtheta(dv)

    # shear srain component
    def gamma_rtheta(du, dv):
        return d_dr(dv) - dv/r + 1.0/r*d_dtheta(du)
    
    epsilon_r(du)
    epsilon_theta(du, dv)
    gamma_rtheta(du, dv)
    

    '''
    S = [[1.0/E, -nu/E, 0.0],
        [-nu/E, 1.0/E, 0.0],
        [0.0, 0.0, 1.0/G]]
    C = df.inv(df.as_matrix(S))
    eps_vector = df.as_vector([epsilon_r(du), epsilon_theta(du, dv), gamma_rtheta(du, dv)])
    sig_vector = dot(C, eps_vector)
    '''

    def sigma_r(du, dv):
        return E/(1.0-nu**2)*( (epsilon_r(du) - alpha*T) + nu*(epsilon_theta(du, dv) - alpha*T) )

    def sigma_theta(du, dv):
        return E/(1.0-nu**2)*( (epsilon_theta(du, dv) - alpha*T) + nu*(epsilon_r(du) - alpha*T) )

    def tau_rtheta(du, dv):
        return G*gamma_rtheta(du, dv)

    # week form
    dF_sigma = - sigma_r(du, dv)*h*r*d_dr(tu)*dx
    dF_sigma = dF_sigma - tau_rtheta(du, dv)*h*d_dtheta(tu)*dx
    dF_sigma = dF_sigma - sigma_theta(du, dv)*h*tu*dx
    dF_sigma = dF_sigma + rho*omega**2*r**2*h*tu*dx

    dF_tau = - tau_rtheta(du, dv)*h*r*d_dr(tv)*dx
    dF_tau = dF_tau - sigma_theta(du, dv)*h*d_dtheta(tv)*dx
    dF_tau = dF_tau + tau_rtheta(du, dv)*h*tv*dx
    dF_tau = dF_tau + rho*omega_velo*r**2*h*tv*dx
    
    dF = dF_sigma + dF_tau

    # residual 
    F = df.action(dF, unks)

    # solve
    df.solve(F == 0, unks, Dbc)

    # splits solution 
    _u, _v = unks.split(True)

    # displacement
    df.File(save_path + 'u.pvd') << _u
    df.File(save_path + 'v.pvd') << _v


    # ======
    # Analyze
    # ======
    # compute stresses
    sigma_r_pro = df.project(sigma_r(_u, _v), V)
    sigma_r_pro.rename('sigma_r [Pa]', 'sigma_r [Pa]')
    df.File(save_path + 'sigma_r.pvd') << sigma_r_pro

    sigma_theta_pro = df.project(sigma_theta(_u, _v), V) 
    sigma_theta_pro.rename('sigma_theta [Pa]', 'sigma_theta [Pa]')
    df.File(save_path + 'sigma_theta.pvd') << sigma_theta_pro

    # compute von Mises stress
    def von_mises_stress(sigma_r, sigma_theta):
        return df.sqrt(sigma_r**2 + sigma_theta**2 - sigma_r*sigma_theta) 

    von_stress_pro = df.project(von_mises_stress(sigma_r(_u, _v), sigma_theta(_u, _v)), V)
    von_stress_pro.rename('von Mises Stress [Pa]', 'von Mises Stress [Pa]')
    df.File(save_path + 'von_mises_stress.pvd') << von_stress_pro

    tau_rtheta_pro = df.project(tau_rtheta(_u, _v), V)
    tau_rtheta_pro.rename('tau_rtheta [Pa]', 'tau_rtheta [Pa]')
    df.File(save_path + 'tau_rtheta.pvd') << tau_rtheta_pro

    # save results to h5
    rfile = df.HDF5File(mesh.mpi_comm(), save_path + 'results.h5', "w")
    rfile.write(_u, "u")
    rfile.write(_v, "v")
    rfile.write(sigma_r_pro, "sigma_r")
    rfile.write(sigma_theta_pro, "sigma_theta")
    rfile.write(von_stress_pro, "von_mises_stress")
    rfile.write(tau_rtheta_pro, "tau_rtheta_pro")
    rfile.close()





