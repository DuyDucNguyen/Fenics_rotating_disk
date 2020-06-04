import dolfin as df
import math

# Implement the model in Mohammad2017
# cartesian coordinate (x, y) 

# Parameters are taken from 
Rint = 0.04 #[m] inner radius
Rext = 0.1 #[m] outer radius
E = 72000*1e6 #[Pa] Young modulus
rho = 2800 #[kg/m^3] density
omega = 15000*2*math.pi/60 #[rad/s] results in 2x MPA
h = 0.005 #[m] disc thickness
nu = 0.3 # Poisson ratio # FAKE

save_path = "main_stress_no_t/"


from create_mesh import create_mesh
from create_mesh_with_holes import create_mesh_with_holes
#mesh = create_mesh(Rint, Rext, 100)
mesh = create_mesh_with_holes(Rint, Rext, Rint/10, 80)


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

df.File(save_path + 'cell_markers.pvd') << facet_markers



# Create mesh and define function space
V = df.FunctionSpace(mesh, "CG", 1)

# Define boundary condition (homogeneous BC)
u0 = df.Constant(0.0)
bc = df.DirichletBC(V, u0, InnerRadius())


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


# Weak form
dres = r*sigma_r(du)*epsilon_r(tu)*df.dx
dres = dres + sigma_theta(du)*tu*df.dx
dres = dres - Fc*tu*df.dx
dres = dres + 1e9*stb_gls


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
















