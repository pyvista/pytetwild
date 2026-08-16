"""Wrapper for fTetWild."""

from .pytetwild import tetrahedralize, tetrahedralize_pv, tetrahedralize_csg  # noqa: F401
from . import _accessor as _accessor

try:
    # Written by setuptools-scm at build time
    from ._version import version as __version__
except ImportError:  # pragma: no cover
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version("pytetwild")
    except PackageNotFoundError:
        __version__ = "unknown"
