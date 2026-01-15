import math
import unittest
import numpy as np
import sys
sys.path.append("..")

from flip import Domain
from flip.core import Beam2D, DofID, CrossSection, Material, Node

class TestBeam2d(unittest.TestCase):
    def setUp(self):
        self.domain = Domain()
        # Nodes
        self.domain.add_node(Node(1, self.domain, coords=[0.0, 0.0, 0.0]))
        self.domain.add_node(Node(2, self.domain, coords=[3.0, 0.0, 0.0]))
        # Material & section
        self.material = Material("mat01", e=210e9)
        self.domain.add_material(self.material)
        self.cs = CrossSection("cs01", a=0.02, iy=8e-6, k=1.0e32)
        self.domain.add_cs(self.cs)
        # Element
        self.element = Beam2D(1, self.domain, nodes=(1, 2), mat="mat01", cs="cs01")
        self.domain.add_element(self.element)
    
    def test_beam2d_getnodedofs(self):
        val = self.element.get_node_dofs(1)
        expected = [DofID.Dx, DofID.Dz, DofID.Ry]
        self.assertEqual(val, expected)
    def test_beam2d_compute_geo(self):
        ans = self.element.compute_geo()
        expected = {'l': 3.0, 'dx': 3.0, 'dz': 0.0}
        self.assertEqual(ans, expected)
    def test_beam2d_has_hinges(self):
        self.assertFalse(self.element.has_hinges())
    def test_beam2d_compute_T(self):
        T = self.element.compute_T()
        self.assertTrue(np.allclose(T, np.eye(6)))
        # test lcs in node 1
        self.domain.nodes[1].update_lcs({'lx': (math.sqrt(3)/2., 0.0, 0.5), 'ly': (0.0, 1.0, 0.0)})
        T = self.element.compute_T()
        ans = T.T @ T
        self.assertTrue(np.allclose(ans, np.eye(6)))
    def test_beam2d_compute_local_stiffness(self):
        kl = self.element.compute_local_stiffness(ret_condensed=False)['answer']
        self.assertEqual(kl.shape, (6, 6))
        self.assertAlmostEqual(kl[0][0], self.material.e * self.cs.a / 3.0)
        self.assertAlmostEqual(kl[0][3], -self.material.e * self.cs.a / 3.0)
        self.assertAlmostEqual(kl[1][1], 12.0 * self.material.e * self.cs.iy / (3.0 ** 3))
        self.assertAlmostEqual(kl[2][2], 4.0 * self.material.e * self.cs.iy / (3.0))
        self.assertAlmostEqual(kl[1][2], -6.0 * self.material.e * self.cs.iy / (3.0 ** 2))
        #test condensed
        self.element.hinges = (False, True)
        kl_condensed = self.element.compute_local_stiffness(ret_condensed=True)['answer']
        self.assertEqual(kl_condensed.shape, (6, 6))
        self.assertAlmostEqual(kl_condensed[0][0], self.material.e * self.cs.a / 3.0)
        self.assertAlmostEqual(kl_condensed[1][1], 3.0 * self.material.e * self.cs.iy / (3.0 ** 3))
        self.assertAlmostEqual(kl_condensed[1][2], -3.0 * self.material.e * self.cs.iy / (3.0 ** 2))
        self.assertAlmostEqual(kl_condensed[5][5], 0.0)

if __name__ == '__main__':
    unittest.main()