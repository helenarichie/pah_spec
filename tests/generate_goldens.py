import os
import sys

import numpy as np

import pah_spec  # pylint: disable=wrong-import-position

test_dir = os.path.dirname(__file__)
golden_path = os.path.join(test_dir, "../data/test_data/pah_spec_golden.csv")
basis_path = os.path.join(test_dir, "../data/test_data/basis_spectra_low_res/")
ps = pah_spec.PahSpec(basis_directory=basis_path)
spectrum_neu, spectrum_ion = ps.generate_spectrum()

golden = np.column_stack((ps.emission_wavelengths.value, spectrum_neu.value, spectrum_ion.value))
# TODO: add command line flag to save goldens to user-specified filename
np.savetxt(
    golden_path,
    golden,
    delimiter=",",
    header="wav [um], p_lambda_neu [erg/(cm s)], p_lambda_ion [erg/(cm s)]",
)
