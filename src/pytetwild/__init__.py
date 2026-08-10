"""Wrapper for fTetWild."""

from . import _accessor as _accessor
from ._version import __version__  # noqa: F401
from .pytetwild import tetrahedralize, tetrahedralize_csg, tetrahedralize_pv  # noqa: F401
