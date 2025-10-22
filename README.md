# `pah_spec`: a code for generating PAH spectra using arbitrary radiation fields

Dependencies: astropy, numpy, pandas, scipy. Tested with python 3.13.7, astropy 6.1.2, numpy 2.3.2, pandas 2.3.2, scipy 1.16.1.

For instructions on code usage, see the examples in [code/generate_spectrum_example.ipynb](https://github.com/helenarichie/pah_spec/blob/main/code/generate_spectrum_example.ipynb) and [code/generate_basis_spectra_example.ipynb](https://github.com/helenarichie/pah_spec/blob/main/code/generate_basis_spectra_example.ipynb).

`pah_spec.PahSpec` requires that a set of basis spectra can be found in the `data/basis_spectra` folder of this repository. We provide a pre-computed set of basis spectra, which can be downloaded from [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LUUXEJ).
