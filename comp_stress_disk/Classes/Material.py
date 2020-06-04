class Material():
    def __init__(self, E, nu, rho):
        self.E = E     # [Pa] Young modulus 
        self.nu = nu   # Poisson ratio
        self.rho = rho # [kg/m^3] density