# `pah_spec`: a code for generating PAH spectra using arbitrary radiation fields

Dependencies: astropy, numpy, pandas, scipy. Tested with python 3.13.7, astropy 6.1.2, numpy 2.3.2, pandas 2.3.2, scipy 1.16.1.

For instructions on code usage, see the examples in [examples/generate_spectrum_example.ipynb](https://github.com/helenarichie/pah_spec/blob/main/examples/generate_spectrum_example.ipynb) and [examples/generate_basis_spectra_example.ipynb](https://github.com/helenarichie/pah_spec/blob/main/examples/generate_basis_spectra_example.ipynb).

`pah_spec.PahSpec` requires that a set of basis spectra can be found in the `data/basis_spectra` folder of this repository. We provide a pre-computed set of basis spectra, which can be downloaded from [here](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/LUUXEJ).

## Installation

Currently, you **MUST** install this package from source. You can skip the rest of this section if you are already familiar with installation of python packages.

There are 3 ways to do this:

1. If you don't plan to download the repository, you can invoke:

   ```sh
   python -m pip install --user git+https://github.com/helenarichie/pah_spec
   ```

   Behind the scenes, pip will download the full repository, install the python package, and delete the repository.

2. If you want to download the repository and then "freeze" the installed package (i.e. the `pah_spec` package won't be affected by any modifications to python files in the repository, unless you perform a fresh installation of the package), you can invoke the following command from the root of the repository:

   ```sh
   python -m pip install --user .
   ```

3. If you want the package to use the most up-to-date versions of the source files in your repository (e.g. for development), you can invoke the following command from the root of the repository:

   ```sh
   python -m pip install --user -e .
   ```

## Development

To run the tests, you should install ``pah_spec`` and then execute the following command from the root of the repository:

```sh
python -m unittest discover -s tests -p "pah_spec_test.py"
```
