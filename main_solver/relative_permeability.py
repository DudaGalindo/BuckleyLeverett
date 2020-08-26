import numpy as np

class RelativePerm:
    """ Class created for the relative permeability calculation """
    def __init__(self, data):
        self.Swr = data['Swr']
        self.Sor = data['Sor']
        self.n_o = data['n_o']
        self.n_w = data['n_w']

    def run(self, Sw):
        krs = np.empty((2,len(Sw)))
        krs[0,:] = self.kro(Sw)
        krs[1,:] = self.krw(Sw)
        return krs

    def kro(self, Sw):
        kro = ((1 - Sw - self.Swr)/(1 - self.Swr - self.Sor))**self.n_o
        return kro

    def krw(self, Sw):
        krw = ((Sw - self.Swr)/(1 - self.Swr - self.Sor))**self.n_w
        return krw
