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


def test_tetrahedralize_matches_direct_api():
    direct = pytetwild.tetrahedralize_pv(pv.Icosphere(nsub=1), edge_length_fac=1.0)
    via_accessor = pv.Icosphere(nsub=1).tetwild.tetrahedralize(edge_length_fac=1.0)
    assert direct.n_cells == via_accessor.n_cells
    assert direct.n_points == via_accessor.n_points


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
