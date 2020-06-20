import comp_stress_disk as csd
from comp_stress_disk.Classes.DiskMesh import DiskMesh
from comp_stress_disk.Classes.Material import Material
from comp_stress_disk.Classes.Simulation import Simulation
from comp_stress_disk.Classes.Analysis import Analysis
from comp_stress_disk.Classes.Sigma_r import Sigma_r
from comp_stress_disk.Classes.Sigma_theta import Sigma_theta
from comp_stress_disk.Classes.Sigma_vm import Sigma_vm
import dolfin as df
import matplotlib.pyplot as plt
import math


# ======
# Mesh Object
# ======
Rint = 0.04 #[m] inner radius
Rext = 0.1 #[m] outer radius
h = 0.005 #[m] disc thickness
Rholes = 0.0
omega = 15000*2*math.pi/60 #[rad/s] results in 2x MPA
res = 80

mesh_obj = DiskMesh(Rint, Rext, Rholes, res)
mesh = mesh_obj.create()




# ======
# Material Object
# ======
# Material properties taken from Mohammad2017
E = 72000*1e6
nu = 0.3 
rho = 2800 #[kg/m^3] density
mat_obj = Material(E, nu, rho)




#TODO: Simulation Object
sim1 = Simulation(mesh_obj = mesh_obj, 
                  mat_obj = mat_obj,
                  axisymmetric = False,
                  heat = False,
                  bc = 'CF',
                  omega = omega, 
                  save_path = 'Simulation_CF/'
                  )
#sim1.run()



#sigma_r = Solution(sim1, 'sigma_r')
sigma_r = Sigma_r(sim1)
#sigma_r.plot()
X, Y = sigma_r.measure_line()
sigma_r.plot_line()

sigma_theta = Sigma_theta(sim1)
sigma_theta.plot_line()

sigma_vm = Sigma_vm(sim1)
sigma_vm.plot_line()

quit()


ana1 = Analysis(sim1)
#ana1.load()
ana1.plot_line(ana1.sigma_r)




quit()

x, sigma_r = ana1.measure_line(ana1.sigma_r)
x, sigma_theta = ana1.measure_line(ana1.sigma_theta)
x, sigma_vm = ana1.measure_line(ana1.sigma_vm)

plt.figure()
plt.plot(x, sigma_r, label="alpha=0")
plt.xlabel('r [m]')
plt.ylabel('sigma_r [Pa]')
plt.legend()

plt.figure()
plt.plot(x, sigma_theta, label="alpha=0")

plt.figure()
plt.plot(x, sigma_vm, label="alpha=0")

plt.show()


quit()


sim2 = Simulation(mesh_obj = mesh_obj, 
                  mat_obj = mat_obj, 
                  bc = 'CC', 
                  omega = omega, 
                  save_path = 'Simulation_CC/')
#sim2.run()


sim3 = Simulation(mesh_obj = mesh_obj, 
                  mat_obj = mat_obj, 
                  bc = 'FF', 
                  omega = omega, 
                  save_path = 'Simulation_FF/')

#sim3.run()


#TODO: Disk Object
