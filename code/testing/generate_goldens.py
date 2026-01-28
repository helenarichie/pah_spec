import astropy.units as u
import numpy as np
import pah_spec # TODO: IMPORT THE LOCAL VERSION!!!!!!

ps = pah_spec.PahSpec()
spectrum_neu, spectrum_ion = ps.generate_spectrum()

golden = np.column_stack((ps.emission_wavelengths.value, spectrum_neu.value, spectrum_ion.value))
# TODO: add command line flag to save goldens to user-specified filename
np.savetxt("../../data/test_data/pah_spec_golden.csv", golden, delimiter=",", 
           header="wav [um], p_lambda,neu [erg/(cm s)], p_lambda,neu [erg/(cm s)]")