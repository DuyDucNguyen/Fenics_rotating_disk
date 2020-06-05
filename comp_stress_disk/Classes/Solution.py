import dolfin as df
import numpy as np
import matplotlib.pyplot as plt


class Solution():
    def __init__(self, simu_obj, name):
        self.simu_obj = simu_obj
        self.name = name
        self.sol = self.load()
    
    def load(self):
        '''load H5 file and read solution'''
        load_path = self.simu_obj.save_path
        mesh = self.simu_obj.mesh_obj.create()
        V = df.FunctionSpace(mesh, 'CG', 1)
        # Load solution
        sol = df.Function(V)
        input_file = df.HDF5File(mesh.mpi_comm(), load_path + 'results.h5', "r")
        input_file.read(sol, self.name)
        input_file.close()
        return sol

    def plot(self):
        '''plot solution'''
        plt.figure()
        fig = df.plot(self.sol)
        plt.xlabel('r [m]')
        plt.ylabel(self.name + ' [Pa]')
        plt.colorbar(fig)
        plt.show()

    def measure_line(self):
        '''measure value from Rint to Rext
        Returns:
            X (list): list of x-coord from Rint to Rext
            Y (list): list of solution(x)
        '''
        Rint = self.simu_obj.mesh_obj.Rint
        Rext = self.simu_obj.mesh_obj.Rext
        X = np.arange(Rint, Rext, (Rext-Rint)/100)
        X = np.append(X, Rext)
        points = [df.Point(x, 0) for x in X]
        Y = [self.sol(p) for p in points]
        return X, Y

    def plot_line(self):
        '''plot line measurement'''
        X, Y = self.measure_line()
        plt.figure()
        fig = plt.plot(X, Y)
        plt.xlabel('r [m]')
        plt.ylabel(self.name + ' [Pa]')
        plt.show()
        return