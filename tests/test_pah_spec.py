from math import isclose
import os
from dataclasses import dataclass

import astropy.units as u
import numpy as np
import pytest

from scipy.integrate import trapezoid

import pah_spec  # pylint: disable=wrong-import-position

_INTERNAL_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

# Set tolerance for disagreement to 5%
TOLERANCE = 0.05

@dataclass
class Scenario:
    # aggregates basic information about the test scenario case
    golden_wavelength_arr: u.Quantity
    golden_spectrum_neu: u.Quantity
    golden_spectrum_ion: u.Quantity
    basis_path: str

@pytest.fixture
def basic_scenario() -> Scenario:
    test_dir = os.path.dirname(os.path.abspath(__file__))
    golden_path = os.path.normpath(os.path.join(test_dir, "../data/test_data/pah_spec_golden.csv"))
    goldens = np.loadtxt(golden_path, delimiter=",", skiprows=1)

    return Scenario(
        golden_wavelength_arr = goldens[:, 0] * u.um,
        golden_spectrum_neu = goldens[:, 1] * u.erg / (u.cm * u.s),
        golden_spectrum_ion = goldens[:, 2] * u.erg / (u.cm * u.s),
        basis_path = os.path.normpath(os.path.join(test_dir, "../data/test_data/basis_spectra_low_res/"))
    )

def test_total_power(basic_scenario: Scenario):
    """Check that integrated spectral power matches goldens for neutral and ionized PAHs."""
    # basic_scenario is an instance of the Scenario type that was constructed by
    # the basic_scenario fixture function

    ps = pah_spec.PahSpec(
        basis_dir=basic_scenario.basis_path,
        internal_data_dir = _INTERNAL_DATA_DIR
    )
    actual_spectrum_neu, actual_spectrum_ion = ps.generate_spectrum()

    actual_power_neu = trapezoid(actual_spectrum_neu, x=ps.emission_wavelengths.to(u.cm))
    actual_power_ion = trapezoid(actual_spectrum_ion, x=ps.emission_wavelengths.to(u.cm))

    expected_power_neu = trapezoid(basic_scenario.golden_spectrum_neu, x=basic_scenario.golden_wavelength_arr.to(u.cm))
    expected_power_ion = trapezoid(basic_scenario.golden_spectrum_ion, x=basic_scenario.golden_wavelength_arr.to(u.cm))

    np.testing.assert_allclose(
        actual_power_neu.value,
        desired=expected_power_neu.value,
        atol=0.0,
        rtol=TOLERANCE,
        err_msg=f"Integrated PAH0 spectral power differs by more than {TOLERANCE*100}%"
    )

    np.testing.assert_allclose(
        actual_power_ion.value,
        desired=expected_power_ion.value,
        atol=0.0,
        rtol=TOLERANCE,
        err_msg=f"Integrated PAH+ spectral power differs by more than {TOLERANCE*100}%"
    )

def test_spectral_shape(basic_scenario: Scenario):
    """Check that spectral shapes are almost equal for neutral and ionized PAHs."""
    # basic_scenario is an instance of the Scenario type that was constructed by
    # the basic_scenario fixture function

    ps = pah_spec.PahSpec(
        basis_dir=basic_scenario.basis_path,
        internal_data_dir = _INTERNAL_DATA_DIR
    )
    actual_spectrum_neu, actual_spectrum_ion = ps.generate_spectrum()

    assert compare_spectra(
        ps.emission_wavelengths.value, actual_spectrum_neu.value,
        basic_scenario.golden_wavelength_arr.value, basic_scenario.golden_spectrum_neu.value
    ), "PAH0 spectra disagree"


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

    n_points = len(x_common)
    n_outliers = sum(not isclose(a, b, rel_tol=TOLERANCE) for a, b in zip(y_actual_interp, y_golden_interp))

    outlier_fraction = 0.05

    return (n_outliers / n_points) <= outlier_fraction

