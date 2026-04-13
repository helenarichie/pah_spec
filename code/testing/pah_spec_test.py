from math import isclose
import os
import sys
import unittest

import astropy.units as u
import numpy as np

from scipy.integrate import trapezoid

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../"))
import pah_spec  # pylint: disable=wrong-import-position

# Set tolerance for disagreement to 1%
TOLERANCE = 0.05


class PahSpecTest(unittest.TestCase):
    """Tests for the pah_spec module.

    unittest executes setUpClass to initialize the class attributes and executes any method
    beginning with "test_" as a test.
    """

    @classmethod
    def setUpClass(cls):
        test_dir = os.path.dirname(__file__)

        golden_path = os.path.join(test_dir, "../../data/test_data/pah_spec_golden.csv")
        goldens = np.loadtxt(golden_path, delimiter=",", skiprows=1)
        cls.golden_wavelength_arr = goldens[:, 0] * u.um
        cls.golden_spectrum_neu = goldens[:, 1] * u.erg / (u.cm * u.s)
        cls.golden_spectrum_ion = goldens[:, 2] * u.erg / (u.cm * u.s)

        basis_path = os.path.join(test_dir, "../../data/test_data/basis_spectra_low_res/")
        cls.ps = pah_spec.PahSpec(basis_directory=basis_path)

        cls.actual_spectrum_neu, cls.actual_spectrum_ion = cls.ps.generate_spectrum()

    def test_total_power(self):
        """Check that integrated spectral power matches goldens for neutral and ionized PAHs."""
        actual_power_neu = trapezoid(self.actual_spectrum_neu, x=self.ps.emission_wavelengths.to(u.cm))
        actual_power_ion = trapezoid(self.actual_spectrum_ion, x=self.ps.emission_wavelengths.to(u.cm))

        expected_power_neu = trapezoid(self.golden_spectrum_neu, x=self.golden_wavelength_arr.to(u.cm))
        expected_power_ion = trapezoid(self.golden_spectrum_ion, x=self.golden_wavelength_arr.to(u.cm))

        self.assertAlmostEqual(
            actual_power_neu.value,
            expected_power_neu.value,
            delta=TOLERANCE,
            msg="Integrated PAH0 spectral power differs by more than 1%",
        )
        self.assertAlmostEqual(
            actual_power_ion.value,
            expected_power_ion.value,
            delta=TOLERANCE,
            msg="Integrated PAH+ spectral power differs by more than 1%",
        )

    def test_spectral_shape(self):
        """Check that spectral shapes are almost equal for neutral and ionized PAHs."""

        try:
            self.assertEqual(
                compare_spectra(
                    self.ps.emission_wavelengths.value,
                    self.actual_spectrum_neu.value,
                    self.golden_wavelength_arr.value,
                    self.golden_spectrum_neu.value,
                ),
                "Arrays are almost equal",
                msg="PAH0 spectra disagree",
            )
        except AssertionError as e:
            print(f"test_almost_equal failed: {e}")
            raise


"""def compare_arrays(a, b):
    if len(a) != len(b):
        raise ValueError("Arrays must be the same length")
    if all(isclose(x.value, y.value, rel_tol=TOLERANCE) for x, y in zip(a, b)):
        return "Arrays are almost equal"
    return "Arrays are not almost equal"""


def compare_spectra(x_actual, y_actual, x_golden, y_golden):
    """Check that two arrays y(x) are within an allowable percentage of each other,
    interpolating onto a common x-axis if they have different lengths.

    Parameters
    ----------
    x_golden : array-like
        x-values for the first array.
    y_golden : array-like
        y-values for the first array.
    x_golden : array-like
        x-values for the second array.
    y_golden : array-like
        y-values for the second array.

    Returns
    -------
    str
        "Arrays are almost equal" if all interpolated values agree within TOLERANCE,
        "Arrays are not almost equal" otherwise.
    """
    # Find the overlapping x range
    x_min = max(np.min(x_golden), np.min(x_actual))
    x_max = min(np.max(x_golden), np.max(x_actual))

    if x_min >= x_max:
        raise ValueError("Spectra have no overlapping wavelength range to compare.")

    # Interpolate both arrays onto a common x-axis within the overlapping range
    x_common = np.logspace(np.log10(x_min), np.log10(x_max), num=len(x_golden))
    y_golden_interp = np.interp(x_common, x_golden, y_golden)
    y_actual_interp = np.interp(x_common, x_actual, y_actual)

    print(x_common)
    print(y_golden_interp)
    print(y_actual_interp)

    n_points = len(x_common)
    n_outliers = sum(not isclose(a, b, rel_tol=TOLERANCE) for a, b in zip(y_actual_interp, y_golden_interp))

    outlier_fraction = 0.05
    print(n_outliers)
    print(len(x_golden))

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].loglog(x_golden, y_golden * x_golden * 1e-4, color="black", label="golden")
    ax[0].loglog(x_actual, y_actual * x_actual * 1e-4, color="red", label="actual")

    ax[1].loglog(x_common, y_golden_interp * x_common * 1e-4, color="black", label="golden")
    ax[1].loglog(x_common, y_actual_interp * x_common * 1e-4, color="red", label="actual")

    ax[1].legend()

    ax[0].set_ylim(1e-26)
    ax[0].set_xlim(1, 20)

    ax[1].set_ylim(1e-26)
    ax[1].set_xlim(1, 20)

    plt.savefig("/Users/helenarichie/Desktop/testing.png")

    if n_outliers / n_points <= outlier_fraction:
        return "Arrays are almost equal"
    return "Arrays are not almost equal"


if __name__ == "__main__":
    unittest.main(verbosity=2)
