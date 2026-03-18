import io
from math import isclose
import os
import unittest
import zipfile

import astropy.units as u
import numpy as np
import requests

from scipy.integrate import trapezoid

import pah_spec

# Set tolerance for disagreement to 1%
TOLERANCE = 0.01


class PahSpecTest(unittest.TestCase):
    """Tests for the pah_spec module.

    unittest executes setUpClass to initialize the class attributes and executes any method
    beginning with "test_" as a test.
    """

    @classmethod
    def setUpClass(cls):

        # Download and save data file before instantiating pah_spec
        url = "https://dataverse.harvard.edu/api/access/dataset/:persistentId/?persistentId=doi:10.7910/DVN/LUUXEJ"
        response = requests.get(url, allow_redirects=True)
        response.raise_for_status()

        extract_path = os.path.join(os.path.dirname(__file__), "../../data/basis_spectra")
        os.makedirs(extract_path, exist_ok=True)

        with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
            zf.extractall(extract_path)

        test_dir = os.path.dirname(__file__)
        golden_path = os.path.join(test_dir, "../../data/test_data/pah_spec_golden.csv")
        goldens = np.loadtxt(golden_path, delimiter=",", skiprows=1)
        cls.golden_wavelength_arr = goldens[:, 0] * u.um
        cls.golden_spectrum_neu = goldens[:, 1] * u.erg / (u.cm * u.s)
        cls.golden_spectrum_ion = goldens[:, 2] * u.erg / (u.cm * u.s)
        cls.ps = pah_spec.PahSpec()

        cls.actual_spectrum_neu, cls.actual_spectrum_ion = cls.ps.generate_spectrum()

        # print(pah_spec.GRAIN_SIZES)

    def test_total_power(self):
        """Check that integrated spectral power matches goldens for neutral and ionized PAHs."""
        actual_power_neu = trapezoid(self.actual_spectrum_neu, x=self.ps.emission_wavelengths.to(u.cm))
        actual_power_ion = trapezoid(self.actual_spectrum_ion, x=self.ps.emission_wavelengths.to(u.cm))

        expected_power_neu = trapezoid(self.golden_spectrum_neu, x=self.golden_wavelength_arr.to(u.cm))
        expected_power_ion = trapezoid(self.golden_spectrum_ion, x=self.golden_wavelength_arr.to(u.cm))

        # TODO assertAlmostEqual
        self.assertAlmostEqual(
            actual_power_neu,
            expected_power_neu,
            delta=TOLERANCE,
            msg="Integrated PAH0 spectral power differs by more than 1%",
        )
        self.assertAlmostEqual(
            actual_power_ion,
            expected_power_ion,
            delta=TOLERANCE,
            msg="Integrated PAH+ spectral power differs by more than 1%",
        )

    def test_spectral_shape(self):
        """Check that spectral shapes are almost equal for neutral and ionized PAHs."""

        try:
            self.assertEqual(
                compare_arrays(self.actual_spectrum_neu, self.golden_spectrum_neu),
                "Arrays are almost equal",
                msg="PAH0 spectra disagree",
            )
        except AssertionError as e:
            print(f"test_almost_equal failed: {e}")
            raise


def compare_arrays(a: list[float], b: list[float]) -> str:
    """Check that each element of two input arrays are within an allowable percentage of each other."""
    if len(a) != len(b):
        raise ValueError("Arrays must be the same length")
    if all(isclose(x.value, y.value, rel_tol=TOLERANCE) for x, y in zip(a, b)):
        return "Arrays are almost equal"
    return "Arrays are not almost equal"


if __name__ == "__main__":
    unittest.main(verbosity=2)
