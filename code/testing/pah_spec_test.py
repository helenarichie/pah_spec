import astropy.units as u
import unittest
import numpy as np
import pah_spec # TODO: IMPORT THE LOCAL VERSION!!!!!!
from scipy.integrate import trapezoid

class PahSpecTest(unittest.TestCase):
    def setUp(self):
        goldens = np.loadtxt('../../data/test_data/pah_spec_golden.csv', delimiter=',', skiprows=1)
        self.golden_wavelength_arr = goldens[:, 0] * u.um
        self.golden_spectrum_neu = goldens[:, 1] * u.erg / (u.cm * u.s)
        self.golden_spectrum_ion = goldens[:, 2] * u.erg / (u.cm * u.s)
        self.ps = pah_spec.PahSpec()
        print(pah_spec.GRAIN_SIZES)

    def test_total_power(self):
        actual_spectrum_neu, actual_spectrum_ion = self.ps.generate_spectrum()

        actual_power_neu = trapezoid(actual_spectrum_neu, x=self.ps.emission_wavelengths.to(u.cm))
        actual_power_ion = trapezoid(actual_spectrum_ion, x=self.ps.emission_wavelengths.to(u.cm))

        expected_power_neu = trapezoid(self.golden_spectrum_neu, x=self.golden_wavelength_arr.to(u.cm))
        expected_power_ion = trapezoid(self.golden_spectrum_ion, x=self.golden_wavelength_arr.to(u.cm))

        # TODO assertAlmostEqual
        self.assertEqual(actual_power_neu, expected_power_neu)
        self.assertEqual(actual_power_ion, expected_power_ion)

if __name__ == '__main__':
    unittest.main()