from .pah_spec import PahSpec, calc_pah_energy, GRAIN_SIZES
from .pah_spec import __doc__

# we consider __version__ to be a public name, but we don't list it the __all__
# variable to avoid importing it when a user writes ``from pah_spec import *``
# -> we import using `as` to silence the ruff warning that __version__ is unused
from ._version import __version__ as __version__

__all__ = ["PahSpec", "calc_pah_energy", "GRAIN_SIZES"]
