"""Tests for the ``.tetwild`` accessor registered on :class:`pyvista.PolyData`."""

from __future__ import annotations

import pytest
import pyvista as pv

import pytetwild  # noqa: F401 — registers the ``.tetwild`` accessor
from pytetwild import _accessor

pytestmark = pytest.mark.skipif(
    not _accessor.HAS_ACCESSOR_REGISTRY,
    reason="requires pyvista >= 0.48 dataset accessor registry",
)


def test_accessor_attached_on_polydata():
    assert hasattr(pv.Sphere(), "tetwild")


def test_accessor_not_attached_on_non_polydata():
    assert not hasattr(pv.ImageData(), "tetwild")
    assert not hasattr(pv.UnstructuredGrid(), "tetwild")


def test_accessor_cached_per_instance():
    sphere = pv.Icosphere(nsub=1)
    assert sphere.tetwild is sphere.tetwild


def test_tetrahedralize_returns_unstructured_grid():
    grid = pv.Icosphere(nsub=1).tetwild.tetrahedralize(edge_length_fac=1.0)
    assert isinstance(grid, pv.UnstructuredGrid)
    assert grid.n_cells > 0
    assert grid.n_points > 0


def test_tetrahedralize_delegates_to_tetrahedralize_pv(monkeypatch):
    """Accessor forwards the parent mesh and kwargs to ``tetrahedralize_pv``.

    fTetWild is non-deterministic across independent runs, so comparing
    cell/point counts of two separate calls is unreliable. The accessor's
    contract is delegation, which is what this asserts.
    """
    sphere = pv.Icosphere(nsub=1)
    sentinel = pv.UnstructuredGrid()
    calls = []

    def fake_tetrahedralize_pv(mesh, **kwargs):
        calls.append((mesh, kwargs))
        return sentinel

    monkeypatch.setattr("pytetwild.pytetwild.tetrahedralize_pv", fake_tetrahedralize_pv)

    result = sphere.tetwild.tetrahedralize(edge_length_fac=1.0, optimize=False)

    assert result is sentinel
    assert len(calls) == 1
    mesh_arg, kwargs_arg = calls[0]
    assert mesh_arg is sphere
    assert kwargs_arg == {"edge_length_fac": 1.0, "optimize": False}


def test_chains_with_core_filters():
    """Accessor result chains cleanly with a core PyVista filter."""
    grid = pv.Icosphere(nsub=1).tetwild.tetrahedralize(edge_length_fac=1.0).extract_cells([0])
    assert isinstance(grid, pv.UnstructuredGrid)
    assert grid.n_cells == 1


def test_registered_record_reports_pytetwild_as_source():
    records = [r for r in pv.registered_accessors() if r.name == "tetwild"]
    assert len(records) == 1
    record = records[0]
    assert record.target is pv.PolyData
    assert record.source.startswith("pytetwild._accessor")
