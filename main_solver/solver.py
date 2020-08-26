import numpy as np
from preprocess.mesh import Mesh
from .relative_permeability import RelativePerm
from .flux import FOUM

class BuckleyLeverett:
    def __init__(self, data):
        self.t = 0.0
        self.Swr = data['Swr']
        self.Sor = data['Sor']
        self.n_el = data['mesh_elements']
        self.L = data['mesh_length']
        self.t_final = data['t_final']
        self.mis = np.empty((2,self.n_el))
        self.mis[0,:] = data['mi_o']
        self.mis[1,:] = data['mi_w']
        self.v = data['v']
        self.krs = RelativePerm(data)
        self.all_results = self.empty_results()

    def __call__(self, data):
        M, Sw = self.initialize(data)
        self.run(M, Sw)
        self.save_results(data)

    def initialize(self, data):
        M = Mesh(self.n_el, data)
        Sw = self.initial_conditions()
        return M, Sw

    def initial_conditions(self):
        Sw = np.zeros(self.n_el)
        Sw[0] = 1 - self.Sor
        Sw[-1] = self.Swr
        return Sw

    def update_mobilities(self, krs):
        ''' Mobilities calculation '''
        lamb = krs/self.mis
        return lamb

    def update_fractional_flux(self, lamb):
        ''' Fractional flow calculations '''
        fs = lamb/np.sum(lamb,axis=1)[:,np.newaxis]
        return fs

    def update_saturation(self, M, Sw, f, delta_t):
        vols_internal = np.array(list(set(M.vols_ID) - set(M.vols_contours_ID)))
        Sw[vols_internal] = Sw[vols_internal] + delta_t * self.v * f.dfw_dSw[0,vols_internal] * \
                            f.dSw_vols[0,vols_internal]
        return Sw

    def update_delta_t(self, f):
        ''' Calculate time step '''
        delta_x = self.L/self.n_el #uniform mesh case
        delta_t = np.min(delta_x/f.dfw_dSw)
        if self.t + delta_t > self.t_final:
            delta_t = self.t_final - self.t
        return delta_t

    def run(self, M, Sw):
        ''' Run simulation '''
        while self.t < self.t_final:
            krs = self.krs.run(Sw)
            lamb = self.update_mobilities(krs)
            fs = self.update_fractional_flux(lamb)
            f = FOUM(M, self.n_el, fs, Sw)
            delta_t = self.update_delta_t(f)
            Sw = self.update_saturation(M, Sw, f, delta_t)
            self.t += delta_t
            self.update_results(Sw)

    def empty_results(self):
        return [np.array(['simulation_time [s]', 'Sw'])]

    def update_results(self, Sw):
        current_results = np.array([self.t, Sw])
        self.all_results.append(current_results)

    def save_results(self, data):
        np.save('results/results_Buckley_Leverett_' + data['mode'] + '_' +
                str(data['mesh_elements']) + '.npy', np.array(self.all_results))
