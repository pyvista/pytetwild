"""Register a ``.tetwild`` accessor on :class:`pyvista.PolyData`.

Importing this module (which :mod:`pytetwild` does on package import)
attaches a :class:`TetWildAccessor` so every :class:`~pyvista.PolyData`
instance exposes the ``.tetwild`` namespace.

Requires PyVista >= 0.48, which introduced ``register_dataset_accessor``.
On older versions importing this module is a no-op.
"""

from __future__ import annotations

from typing import Any

import pyvista as pv


HAS_ACCESSOR_REGISTRY = hasattr(pv, "register_dataset_accessor")


def _register(cls):
    """Register ``cls`` as the ``.tetwild`` accessor when supported.

    Parameters
    ----------
    cls : type
        Accessor class to register.

    Returns
    -------
    type
        ``cls``, registered when PyVista exposes the accessor registry.
    """
    if HAS_ACCESSOR_REGISTRY:
        return pv.register_dataset_accessor("tetwild", pv.PolyData)(cls)
    return cls


@_register
class TetWildAccessor:
    """Tetrahedralization accessor backed by fTetWild.

    Wraps :func:`pytetwild.tetrahedralize_pv` so a user can tetrahedralize
    a surface mesh in a single chained call::

        import pyvista as pv
        import pytetwild  # noqa: F401 — registers the ``.tetwild`` accessor

        grid = pv.Icosphere(nsub=2).tetwild.tetrahedralize(edge_length_fac=0.1)

    The accessor is scoped to :class:`~pyvista.PolyData` because fTetWild
    consumes a triangulated surface. Non-triangular input is
    auto-triangulated by the underlying function with a warning.

    Parameters
    ----------
    mesh : pyvista.PolyData
        Surface mesh exposing the ``.tetwild`` namespace.
    """

    def __init__(self, mesh: pv.PolyData) -> None:
        """Bind the accessor to its parent :class:`~pyvista.PolyData`.

        Parameters
        ----------
        mesh : pyvista.PolyData
            Surface mesh exposing the ``.tetwild`` namespace.
        """
        self._mesh = mesh

    def tetrahedralize(self, **kwargs: Any) -> pv.UnstructuredGrid:
        """Tetrahedralize the surface mesh and return the volume grid.

        Thin wrapper around :func:`pytetwild.tetrahedralize_pv`.
        Forwards all keyword arguments; see that function's docstring
        for the full parameter list (``edge_length_fac``, ``optimize``,
        ``epsilon``, ``stop_energy``, and so on).

        Parameters
        ----------
        **kwargs
            Forwarded to :func:`pytetwild.tetrahedralize_pv`.

        Returns
        -------
        pyvista.UnstructuredGrid
            Tetrahedral volume mesh.
        """
        from pytetwild.pytetwild import tetrahedralize_pv

        return tetrahedralize_pv(self._mesh, **kwargs)
