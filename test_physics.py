import unittest
import physics


class TestPhysics(unittest.TestCase):
    def test_constants(self):
        self.assertEquals(g, 9.81)
        self.assertEquals(density_water, 1000)

    def test_calculate_buoyancy(self):
        self.assertRaises(ValueError, physics.calculate_buoyancy(-5, 6))
        self.assertRaises(ValueError, physics.calculate_buoyancy(6, -5))
        self.assertEquals(physics.calculate_buoyancy(1, density_water), 9810)
        self.assertEquals(physics.calculate_buoyancy(60, density_water), 588600)

    def test_will_it_float(self):
        self.assertRaises(
            ValueError,
            physics.will_it_float(-5, 5),
        )
        self.assertRaises(ValueError, physics.will_it_float(5, -5))
        self.assertTrue(physics.will_it_float(50, 5))
        self.assertFalse(physics.will_it_float(9, 10000))
        self.assertFalse(
            physics.will_it_float(10, 10000)
        )  # = density water, so neither floats nor sinks

    def test_calculate_pressure(self):
        self.assertRaises(ValueError, physics.calculate_pressure(-10))
        self.assertRaises(ValueError, physics.calculate_pressure(0))
        self.assertEquals(physics.calculate_pressure(1), 9810)
        self.assertEquals(physics.calculate_pressure(3657), 35875170)


if __name__ == "__main__":
    unittest.main()
