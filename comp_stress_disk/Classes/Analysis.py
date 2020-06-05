import dolfin as df
import numpy as np
import matplotlib.pyplot as plt





simu_obj = -1

class Analysis():
    def load(self):
        '''
        load H5 file and read u, sigma_r, sigma_theta
        '''
        load_path = self.simu_obj.save_path
        mesh = self.simu_obj.mesh_obj.create()
        V = df.FunctionSpace(mesh, 'CG', 1)
        # Load solution
        u = df.Function(V)
        sigma_r = df.Function(V)
        sigma_theta = df.Function(V)
        vm_stress = df.Function(V)
        input_file = df.HDF5File(mesh.mpi_comm(), load_path + 'results.h5', "r")
        input_file.read(u, "u")
        input_file.read(sigma_r, "sigma_r")
        input_file.read(sigma_theta, "sigma_theta")
        input_file.read(vm_stress, "von_mises_stress")
        input_file.close()
        return u, sigma_r, sigma_theta, vm_stress

    def __init__(self, 
                 simu_obj = simu_obj, 
                ):
        self.simu_obj = simu_obj
        u, sigma_r, sigma_theta, vm_stress = self.load()
        self.u = u
        self.sigma_r = sigma_r
        self.sigma_theta = sigma_theta
        self.sigma_vm = vm_stress

    def measure_line(self, f):
        '''
        measure value from Rint to Rext
        '''
        Rint = self.simu_obj.mesh_obj.Rint
        Rext = self.simu_obj.mesh_obj.Rext
        X = np.arange(Rint, Rext, (Rext-Rint)/100)
        # TODO: add Rext
        #X = X.add(Rext)
        points = [df.Point(x, 0) for x in X]
        measure = [f(p) for p in points]
        return X, measure

    def plot_line(self, f):
        X, measure = self.measure_line(f)
        plt.figure()
        plt.plot(X, measure)
        plt.xlabel('r [m]')
        if f == self.sigma_r:
            plt.ylabel('sigma_r [Pa]')
        elif f == self.sigma_theta: 
            plt.ylabel('sigma_theta [Pa]')
        elif f == self.sigma_vm: 
            plt.ylabel('sigma_vm [Pa]')
        #plt.legend()
        plt.show()
        return