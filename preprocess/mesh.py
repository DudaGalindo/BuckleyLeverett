import numpy as np

class Mesh:
    """ 1-D mesh creation """

    def __init__(self, n_el, data):
        self.vols_ID = np.arange(n_el)
        self.vols_coord = np.linspace(0+1/n_el,1-data['mesh_length']/n_el,n_el)
        self.internal_faces_ID = np.arange(n_el-1)
        self.vols_contours_ID = np.array([0,n_el-1])
        self.neighboring_vols(n_el)

    def neighboring_vols(self, n_el):
        self.conec = np.empty((n_el-1,2),dtype=int)
        self.conec[:,0] = np.arange(n_el-1)
        self.conec[:,1] = np.arange(1,n_el)
