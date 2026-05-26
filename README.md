# pah_spec

``pah_spec`` is a tool for fast, flexible computation of PAH emission spectra with the single photon approximation, as described in [Richie & Hensley (2026)](https://ui.adsabs.harvard.edu/abs/2025arXiv251016861R/abstract).

## Getting ``pah_spec``

Installing ``pah_spec`` is as simple as invoking the following pair of commands:

```sh
# install the pah_spec package
$ python -m pip install --user git+https://github.com/helenarichie/pah_spec

# install the required data files. This can take some time. If the tqdm package
# is installed, progressbars are shown. You can also skip this step and install
# these files programmatically
$ python -m pah_spec download --to-cache all
```

Out installation guide provides more detailed instructions (and describes other installation-related choices).

## Examples:

For instructions on code usage, see the examples in [examples/generate_spectrum_example.ipynb](https://github.com/helenarichie/pah_spec/blob/main/examples/generate_spectrum_example.ipynb) and [examples/generate_basis_spectra_example.ipynb](https://github.com/helenarichie/pah_spec/blob/main/examples/generate_basis_spectra_example.ipynb).
