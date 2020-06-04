import dolfin as df
from ufl import nabla_div
import sys
import os
import numpy as np
import math

from fenics_convert_unit import convert_unit


# Implement the model for ratating annular disc
# Mohammad2017 - Elastic Stress Analysis of Rotating Functionally Graded Annular Disk of Variable Thickness Using Finite Difference Method
# the model is rewritten to fix cartesian coordinate (x, y)


def material_function(target_mesh, cells_list, coeffs):
    """
    turn a coefficient that varies by different materials into a function
    """
    coeff_func = df.Function(df.FunctionSpace(target_mesh, "DG", 0))
    markers = np.asarray(cells_list.array(), dtype=np.int32)
    coeff_func.vector()[:] = np.choose(markers, coeffs)
    return coeff_func


# =========================== Terminal variables ===========================
load_path = 'Geometry/'
mesh_name = 'Circle2d_mohamed.xdmf'
mesh_name = 'Circle2d_mohamed_with_holes.msh'
mesh_name = 'Circle2d_mohamed_with_holes.xdmf'
save_path = 'Circle2d_mohamed_with_holes/'


E_lam = 70000  # lamination Young's modulus
nu_lam = 0.28  # lamination Poisson ratio
rho_lam = 2800  # lamination mass
E_mag = 70000  # magnetic Young's modulus
nu_mag = 0.3  # magnetic Poisson ratio
rho_mag = 2800  # magnetic mass
omega = 15000*2*math.pi  # [Hz] frequency


Rint = 0.04  # [m] inner radius
Rext = 0.1  # [m] outer radius
h0 = 0.005  # [m] disc thickness

#alpha = -2 # match with alpha = -1 in Mohammad2017
#alpha = -1 # match with alpha = 0
#alpha = 0 # match with alpha = 1
alpha = 0
alpha = 1

'''
print('E_lam', E_lam)
print('E_mag', E_mag)
print('nu_lam', nu_lam)
print('nu_mag', nu_mag)
print('rho_lam', rho_lam)
print('rho_mag', rho_mag)
'''

#print('Rint', Rint)
#print('Rext', Rext)


# =========================== Change Units ===========================
# Lame coefficient for constitutive relation
def mu_func(E, nu):
    return E / 2.0 / (1.0 + nu)


def lmbda_func(E, nu):
    return E * nu / (1.0 + nu) / (1.0 - 2.0 * nu)


def E_func(mu, lmbda):
    return mu * (3.0 * lmbda + 2.0 * mu) / (lmbda + mu)


def nu_func(mu, lmbda):
    return lmbda / (2.0 * (lmbda + mu))


## ==== convert unit ====

s_unit = "m"
t_unit = "s"
w_unit = "kg"


E_lam = (
    E_lam
    * convert_unit(1.0, "kg", w_unit)
    / (convert_unit(1.0, "m", s_unit) * convert_unit(1.0, "s", t_unit) ** 2)
)
rho_lam = (
    rho_lam * convert_unit(1.0, "kg", w_unit) / (convert_unit(1.0, "m", s_unit) ** 3)
)


E_mag = (
    E_mag
    * convert_unit(1.0, "kg", w_unit)
    / (convert_unit(1.0, "m", s_unit) * convert_unit(1.0, "s", t_unit) ** 2)
)
rho_mag = (
    rho_mag * convert_unit(1.0, "kg", w_unit) / (convert_unit(1.0, "m", s_unit) ** 3)
)


# =========================== Set Paths ===========================
head, tail = mesh_name.split(".")
# save_path = os.path.splitext(__file__)[0] + '/'
# check_dir(save_path)


# =========================== Load Mesh ===========================
from fenics_meshio_convert_to_xdmf import meshio_convert_to_xdmf

meshio_convert_to_xdmf(load_path, mesh_name, load_path)


# =========================== Load Mesh ===========================
from fenics_meshio_read_xdmf import meshio_read_xdmf

mesh_name = head + ".xdmf"

# read mesh
mesh, cell_markers, facet_markers, tag_map = meshio_read_xdmf(load_path, mesh_name)
is_facet_markers_manual = False

#print(tag_map)

if bool(tag_map): 
    lam_tag = tag_map['annular'][0]
    mag_tag = tag_map['annular'][0]
else: 
    lam_tag = 3
    mag_tag = 3
#print(lam_tag)
#print(mag_tag)


# dimention
dim = mesh.topology().dim()


"""
from mshr import *
# Create circles as Circle(Center, Radius)
outer_circle = Circle(df.Point(0,0), b)
inner_circle = Circle(df.Point(0,0), a)
domain = outer_circle - inner_circle
r = 100   # Resolution of mesh
mesh = generate_mesh(domain, r)
x, y = df.SpatialCoordinate(mesh)
"""

# Define coordinate variable
x = df.Expression("x[0]", degree=1)
y = df.Expression("x[1]", degree=1)


# =========================== Manually mark facet_marker ===========================
class OuterRadius(df.SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0]**2 + x[1]**2 >= Rext**2 - 1e-8

class InnerRadius(df.SubDomain):
    def inside(self, x, on_boundary):
        return x[0]**2 + x[1]**2 < Rint**2 + 1e-8 and on_boundary


if facet_markers == False:
    facet_markers = df.MeshFunction('size_t', mesh, dim-1)
    facet_markers.set_all(0)
    InnerRadius().mark(facet_markers, 1)
    is_facet_markers_manual = True


'''
# =========================== Fix Cell_markers ===========================
# WARNING: only use for new_machine.gmsh
for c in df.cells(mesh):
    if cell_markers[c] > 1: 
        cell_markers[c] = 1
'''


# save marker functions into files readable with Paraview
df.File(save_path + "Marker_Functions/" + "cell_markers.pvd") << cell_markers
df.File(save_path + "Marker_Functions/" + "facet_markers.pvd") << facet_markers


# =========================== Define measurement ===========================
# define subdomains and cell measurement
dx = df.Measure("dx", domain=mesh, subdomain_data=cell_markers)

# Measurement associated with the Exterior boundaries
ds = df.Measure("ds", domain=mesh, subdomain_data=facet_markers)

# define interface and facet measurement
dS = df.Measure("dS", domain=mesh, subdomain_data=facet_markers)


# =========================== Define Material Expression ===========================
class SET_PARAMETER(df.UserExpression):
    def __init__(self, subdomains, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == lam_tag:  # Lamination tag
            values[0] = self.k_0
        elif self.subdomains[cell.index] == mag_tag: # magnet tag
            values[0] = self.k_1

    def value_shape(self):
        return ()


E = SET_PARAMETER(subdomains=cell_markers, k_0=E_lam, k_1=E_mag, degree=0)
E_pro = df.project(E, df.FunctionSpace(mesh, "DG", 0))
E_pro.rename("E", "E")
df.File(save_path + "Material/" + "E.pvd") << E_pro


rho = SET_PARAMETER(subdomains=cell_markers, k_0=rho_lam, k_1=rho_mag, degree=0)
rho_pro = df.project(rho, df.FunctionSpace(mesh, "DG", 0))
rho_pro.rename("rho", "rho")
df.File(save_path + "Material/" + "rho.pvd") << rho_pro


nu = SET_PARAMETER(subdomains=cell_markers, k_0=nu_lam, k_1=nu_mag, degree=0)
nu_pro = df.project(nu, df.FunctionSpace(mesh, "DG", 0))
nu_pro.rename("nu", "nu")
df.File(save_path + "Material/" + "nu.pvd") << nu_pro


mu = mu_func(E, nu)
lmbda = lmbda_func(E, nu)
m_c = rho


# =========================== Define Finite Element ===========================
V = df.FunctionSpace(mesh, "CG", 1)

# Define variational problem
du = df.TrialFunction(V)
tu = df.TestFunction(V)
# displacement in radial direction u(x,y)
u = df.Function(V, name="displacement")


# =========================== Dirichlet Boundary Conditions ===========================
#if not is_facet_markers_manual: 
#    bc = df.DirichletBC(V, df.Constant(0), facet_markers, tag_map["interior"][0])
#else: 
#    bc = df.DirichletBC(V, df.Constant(0), facet_markers, 1)

bc = df.DirichletBC(V, df.Constant(0), facet_markers, 1)


# =========================== Define the model ===========================
class THETA(df.UserExpression):
    def eval(self, values, x):
        values[0] = math.atan2(x[1], x[0])

    def value_shape(self):
        # return (1,) # vector
        return ()  # scalar


theta = THETA(degree=1)
# theta_int = df.interpolate(theta, df.FunctionSpace(mesh, "DG", 0))
# df.File(save_path + 'theta.pvd') << theta_int


class RADIUS(df.UserExpression):
    def eval(self, values, x):
        values[0] = df.sqrt(x[0] * x[0] + x[1] * x[1])

    def value_shape(self):
        return ()  # scalar


#r = df.sqrt(x**2 + y**2)
r = RADIUS(degree=1)


# strain radial
def epsilon_r(du):
    return 1.0 / r * (x * df.Dx(du, 0) + y * df.Dx(du, 1))

    
# strain circumferential
def epsilon_theta(du):
    return du / df.sqrt(x ** 2 + y ** 2)


# radial stress # train-stress relation
def sigma_r(du):
    return E / (1.0 - nu ** 2) * (epsilon_r(du) + nu * epsilon_theta(du))


# circumferential stress
def sigma_theta(du):
    return E / (1.0 - nu ** 2) * (nu * epsilon_r(du) + epsilon_theta(du))


# thickness 
t = h0*(r/Rext)**alpha

# define plane stress equilibrium equation
Fc = t*rho*omega**2*r**2


# ================= Streamline Stab: Galerkin Least Square (GLS) Stabilization =================
# stabilization term 

# Cellsize function
cs_func = 2.0*df.Circumradius(mesh) 
cs_pro = df.project(cs_func, V)
df.File(save_path + 'cell_size.pvd') << cs_pro


# scale of cell size
cs_scale = 1e-3
cs = cs_scale*cs_func


diff_coef = - E*r/(-nu**2+1)
conv_coef = - E/(-nu**2*r+1)
reac_coef = E/(-nu**2*r+r)

tau = pow(2*conv_coef/cs + 4*diff_coef/cs**2 + reac_coef, -1)


def d_dr(du): 
    return 1.0/r*(x*df.Dx(du, 0) + y*df.Dx(du, 1))

def d_dtheta(du):
    return -y*df.Dx(du, 0) + x*df.Dx(du, 1)

def d_dx(du): 
    return df.Dx(du, 0)

def d_dy(du): 
    return df.Dx(du, 1)


def d2_dr2(du): 
    d2u_dx2 = df.Dx(df.Dx(du, 0), 0) 
    du_dxdy = df.Dx(df.Dx(du, 0), 1) 
    d2u_dy2 = df.Dx(df.Dx(du, 1), 1)
    du_dx = d_dx(du)
    du_dy = d_dy(du)
    d2u_dr2 = (x*(-x*(x*du_dx + y*du_dy)/(x**2 + y**2)**(3/2) + (x*d2u_dx2 + y*du_dxdy + du_dx)/df.sqrt(x**2 + y**2)) + y*(-y*(x*du_dx + y*du_dy)/(x**2 + y**2)**(3/2) + (x*du_dxdy + y*d2u_dy2 + du_dy)/df.sqrt(x**2 + y**2)))/df.sqrt(x**2 + y**2)
    return d2u_dr2

def d2_dtheta2(du):
    d2u_dy2 = d_dy(d_dy(du))
    du_dxdy = d_dy(d_dx(du))
    d2u_dx2 = d_dx(d_dx(du))
    return x*(x*d2u_dy2 - y*du_dxdy - d_dx(du)) - y*(x*du_dxdy - y*d2u_dx2 + d_dy(du))



def polar_div(du): 
    return d2_dr2(du) + 1/r*d_dr(du) + 1/r**2*d2_dtheta2(du)


def L(du): 
    return diff_coef*d2_dr2(du) + conv_coef*d_dr(du) + reac_coef*du  + polar_div(du) #+ diff_coef*d2_dtheta2(du) + 
    #return diff_coef*polar_div(du) + conv_coef*d_dr(du) + conv_coef*d_dtheta(du) + reac_coef*du 
    #return polar_div(du) + d_dr(du) + d_dtheta(du) + du #+ diff_coef*d2_dtheta2(du) + 



def Res(du):
    return L(du) - rho*omega**2*r**2



# Galerkin Least Square stabilization term
stb_gls = L(tu)*tau*Res(du)*df.dx 





# Galerkin Least Square stabilization term
stb_gls = L(tu)*tau*Res(du)*df.dx 
#stb_gls = L(u)*tu*df.dx #
#stb_gls = L(tu)*tau*Res(u)*df.dx 


# Weak form
dres = - t*r*sigma_r(du)*epsilon_r(tu)*df.dx
dres = dres - t*sigma_theta(du)*tu*df.dx
dres = dres + Fc*tu*df.dx
dres = dres - stb_gls

# residual
res = df.action(dres, u) 



# Compute solution
df.solve(res == 0, u, bc)


# displacement
df.File(save_path + "displacement.pvd") << u


# define a line
line_x = np.arange(Rint, Rext, (Rext-Rint)/100)
line_p = [df.Point(x, 0.0) for x in line_x]

# compute sigma_r
sigma_r_pro = df.project(sigma_r(u), V)
sigma_r_pro.rename("sigma_r [Pa]", "sigma_r [Pa]")
df.File(save_path + "sigma_r.pvd") << sigma_r_pro
out_val = np.asarray([(x, sigma_r_pro(df.Point(x, 0.0))) for x in line_x])
filename = 'sigma_r_alpha_{}.txt'.format(str(alpha))
np.savetxt(save_path + filename, out_val, delimiter='\t')

# compute sigma_theta
sigma_theta_pro = df.project(sigma_theta(u), V)
sigma_theta_pro.rename("sigma_theta [Pa]", "sigma_theta [Pa]")
df.File(save_path + "sigma_theta.pvd") << sigma_theta_pro
out_val = np.asarray([(x, sigma_theta_pro(df.Point(x, 0.0))) for x in line_x])
filename = 'sigma_theta_alpha_{}.txt'.format(str(alpha))
np.savetxt(save_path + filename, out_val, delimiter='\t')


# compute von Mises stress
def von_mises_stress(sigma_r, sigma_theta):
    return df.sqrt(sigma_r ** 2 + sigma_theta ** 2 - sigma_r * sigma_theta)

von_stress_pro = df.project(von_mises_stress(sigma_r(u), sigma_theta(u)), V)
von_stress_pro.rename("von Mises Stress [Pa]", "von Mises Stress [Pa]")
df.File(save_path + "von_mises_stress.pvd") << von_stress_pro


