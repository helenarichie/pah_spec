"""
pah_spec is a package for generating PAH emission spectra for arbitrary input radiation fields using the single photon
approximation, as described in Richie & Hensley (2026).

The primary functionality of pah_spec lies in the :py:class:`PahSpec` class, which uses a pre-existing set of basis spectra to
execute the :py:func:`PahSpec.generate_spectrum` method. This method takes a radiation field and dust grain size distributions as
(optional) inputs and returns an emission spectrum for neutral and ionized PAHs.

In addition to generating emission spectra, pah_spec contains useful functions for the PAH energy and absorption
cross-sections (:py:func:`calc_pah_energy` and :py:meth:`PahSpec.calc_c_abs`, respectively) using the PAH model from Draine & Li (2001)
and Draine et al. (2021).
"""

from ._core import PahSpec, calc_pah_energy

# we consider __version__ to be a public name, but we don't list it the __all__
# variable to avoid importing it when a user writes ``from pah_spec import *``
# -> we import using `as` to silence the ruff warning that __version__ is unused
from ._version import __version__ as __version__

__all__ = ["PahSpec", "calc_pah_energy"]
