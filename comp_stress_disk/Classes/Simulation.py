from ..Methods.Simulation.comp_axisymmetric_dirichlet import comp_axisymmetric_dirichlet
from ..Methods.Simulation.comp_axisymmetric_pure_neumann import comp_axisymmetric_pure_neumann
from ..Methods.Simulation.comp_full_model_heat_dirichlet import comp_full_model_heat_dirichlet


mesh_obj = -1
mat_obj = -1


def is_bool(val):
    try:
        return {"true":True,"false":False}[str(val).lower()]
    except KeyError:
        print("Invalid input please enter True or False!")


class Simulation():
    def __init__(self, 
                mesh_obj = mesh_obj, 
                mat_obj = mat_obj, 
                axisymmetric = True,
                heat = None,
                bc = 'CF', 
                omega = 0.0,
                save_path = ''
                ):
        
        is_bool(axisymmetric)
        
        self.mat_obj = mat_obj
        self.mesh_obj = mesh_obj
        self.axisymmetric = axisymmetric
        self.heat = heat
        self.bc = bc
        self.omega = omega
        self.save_path = save_path

    #TODO: model: axisymmetry, full, heat, heat axisymmetry
    
    def run(self):
        if self.axisymmetric is True:
            if self.heat is True: 
                pass
            elif self.heat is None:
                if self.bc == 'CF' or self.bc == 'CC':
                    comp_axisymmetric_dirichlet(self.mat_obj, self.mesh_obj, self.bc, self.omega, self.save_path)
                elif self.bc == 'FF':
                    comp_axisymmetric_pure_neumann(self.mat_obj, self.mesh_obj, self.bc, self.omega, self.save_path)
        elif self.axisymmetric is False: 
            comp_full_model_heat_dirichlet(self.mat_obj, self.mesh_obj, self.bc, self.omega, self.save_path)
        return