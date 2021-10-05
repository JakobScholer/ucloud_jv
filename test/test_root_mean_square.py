import unittest
from src.root_mean_square import root_mean_square


class Test(unittest.TestCase):
    # Normal behaviour test case
    def test_root_mean_square_normal(self):
        energy_curve1 = [0.0, 6.790366, 27.059802, 58.906111, 100.867318, 149.478833, 122.747973, 114.464051]
        energy_curve2 = [0.0, 7.023074, 27.118542, 57.480826, 95.409624, 137.375519, 152.960128, 121.460638, 113.663662, 112.157104]
        actual = root_mean_square(energy_curve1, energy_curve2)
        expected = 2.9530847717119375
        self.assertEqual(actual, expected)

    # Two identical curves
    def test_root_mean_square_same(self):
        energy_curve1 = [0.0, 6.790366, 27.059802, 58.906111, 100.867318, 149.478833, 122.747973, 114.464051]
        actual = root_mean_square(energy_curve1, energy_curve1)
        expected = 0
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
