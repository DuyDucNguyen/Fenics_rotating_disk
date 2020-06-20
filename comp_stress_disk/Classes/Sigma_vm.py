from .Solution import Solution

class Sigma_vm(Solution):
    def __init__(self, simu_obj):
        Solution.__init__(self, simu_obj, 'von_mises_stress')

