from ..Methods.Simulation.comp_axisymmetric_dirichlet import comp_axisymmetric_dirichlet
from ..Methods.Simulation.comp_axisymmetric_pure_neumann import comp_axisymmetric_pure_neumann

mesh_obj = -1
mat_obj = -1

class Simulation():
    def __init__(self, 
                mesh_obj = mesh_obj, 
                mat_obj = mat_obj, 
                bc = 'CF', 
                omega = 0.0,
                save_path = ''
                ):
        self.mat_obj = mat_obj
        self.mesh_obj = mesh_obj
        self.bc = bc
        self.omega = omega
        self.save_path = save_path

    #TODO: model: axisymmetry, full, heat, heat axisymmetry
    def run(self):
        if self.bc == 'CF' or self.bc == 'CC':
            comp_axisymmetric_dirichlet(self.mat_obj, self.mesh_obj, self.bc, self.omega, self.save_path)
        elif self.bc == 'FF':
            comp_axisymmetric_pure_neumann(self.mat_obj, self.mesh_obj, self.bc, self.omega, self.save_path)
        return