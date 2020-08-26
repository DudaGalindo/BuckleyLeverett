from main_solver.solver import BuckleyLeverett
import unittest

class Test_BL(unittest.TestCase):
    def test_8_el(self):
        data = dict()
        # INPUT DATA FILE
        data['mesh_elements'] = 8
        data['mesh_length'] = 1
        data['mode'] = 'FOUM'

        # RELATIVE PERMEABILITY DATA
        data['Sor'] = 0
        data['Swr'] = 0
        data['Sw'] = 0
        data['n_o'] = 2
        data['n_w'] = 2

        # VISCOSITY DATA
        data['mi_o'] = 1
        data['mi_w'] = 1

        # TOTAL VELOCITY
        data['v'] = 1

        # TOTAL SIMULATION TIME
        data['t_final'] = 10
        BL = BuckleyLeverett(data)
        BL(data)
