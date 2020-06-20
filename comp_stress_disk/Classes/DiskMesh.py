import dolfin as df
import matplotlib.pyplot as plt
from ..Methods.DiskMesh.create_simple_mesh import create_simple_mesh
from ..Methods.DiskMesh.create_mesh_with_holes import create_mesh_with_holes





class DiskMesh:
    """
        mesh object
        :param Rint: interior raidus
        :param Rext: exterior raidus
        :param res: mesh resolution
        :param Rholes: holes radius
        
    """
    def __init__(self, Rint, Rext, Rholes, res):
        self.Rint = Rint
        self.Rext = Rext
        self.res = res
        self.Rholes = Rholes

    def create(self):
        if self.Rholes == 0:
            mesh = create_simple_mesh(self.Rint, self.Rext, self.res)
        if self.Rholes > 0:
            mesh = create_mesh_with_holes(self.Rint, self.Rext, self.Rholes, self.res)
        return mesh
    

    def save(self, save_path, mesh_name):
        df.File(save_path + mesh_name) << self.create()
        return

    def plot(self):
        plt.figure()
        df.plot(self.create())
        plt.show()

