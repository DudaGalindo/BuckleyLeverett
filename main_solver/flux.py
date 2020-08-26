import numpy as np
import scipy.sparse as sp

class FOUM:
    """ Class created for calculating flux at the interface using First Order Upwind Method """

    def __init__(self, M, n_el, fs, Sw):
        fs_face, Sw_face = self.upwind(M, fs, Sw)
        self.dfw_dSw, self.dSw_vols = self.flux_vols(M, n_el, fs_face, Sw_face)

    def upwind(self, M, fs, Sw):
        """ Approximating fractional flux and water saturation at the interface according to
        the flux direction """
        fs_face_neig = fs[:,M.conec]
        Sw_face_neig = Sw[M.conec]
        fs_face = np.empty_like(fs[:,M.conec[:,0]])
        Sw_face = np.empty_like(Sw_face_neig[:,0])

        """ Approximating phase fractional flux at block faces, fs[0] - oil phase e fs[1] - water phase """
        fs_face[fs_face_neig[:,:,1]>=fs_face_neig[:,:,0]] = fs_face_neig[
                                    fs_face_neig[:,:,1]>=fs_face_neig[:,:,0],1]
        fs_face[fs_face_neig[:,:,1]<fs_face_neig[:,:,0]] = fs_face_neig[
                                    fs_face_neig[:,:,1]<fs_face_neig[:,:,0],0]

        """ Approximating water saturation at block faces """
        Sw_face[fs_face_neig[1,:,1]>=fs_face_neig[1,:,0]] = Sw_face_neig[
                                    fs_face_neig[1,:,1]>=fs_face_neig[1,:,0],1]
        Sw_face[fs_face_neig[1,:,1]<fs_face_neig[1,:,0]] = Sw_face_neig[
                                    fs_face_neig[1,:,1]<fs_face_neig[1,:,0],0]
        return fs_face, Sw_face

    def flux_vols(self, M, n_el, fs_face, Sw_face):
        """ Calculate dfw/dSw and dSw at each grid block """

        ''' fractional flux variation at each control volume '''
        cph = np.arange(2) #2 phases
        lines = np.array([np.repeat(cph,len(M.conec[:,0])), np.repeat(cph,
                                len(M.conec[:,1]))]).astype(int).flatten()
        cols = np.array([np.tile(M.conec[:,0],2), np.tile(M.conec[:,1], 2)]).flatten()
        data = np.array([-fs_face, fs_face]).flatten()
        dfs_vols = sp.csc_matrix((data, (lines, cols)), shape = (2, n_el)).toarray()

        ''' Saturation variation at each control volume '''
        lines = np.array([np.repeat(0,len(M.conec[:,0])), np.repeat(0,
                                    len(M.conec[:,1]))]).astype(int).flatten()
        cols = np.array([M.conec[:,0], M.conec[:,1]]).flatten()
        data = np.array([-Sw_face, Sw_face]).flatten()
        dSw_vols = sp.csc_matrix((data, (lines, cols)), shape = (1,n_el)).toarray()

        np.seterr(all = 'ignore', divide = 'ignore')
        dfw_dSw = dfs_vols[1,:] / dSw_vols
        dfw_dSw[abs(dSw_vols)<1e-15] = 0
        return dfw_dSw, dSw_vols

class MUSCL:
    pass
