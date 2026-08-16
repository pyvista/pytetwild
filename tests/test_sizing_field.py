"""Tests for the input sizing field (fTetWild's ``--bg-mesh``)."""

from __future__ import annotations

import numpy as np
import pyvista as pv
import pytest

from pytetwild import tetrahedralize, tetrahedralize_pv

# Sphere of diameter 1, so a bounding box diagonal of about 1.73 and a
# default ideal edge length of about 0.087.
FINE = 0.04
COARSE = 0.2

# Ratio of mean cell volume between the coarse and fine halves. Over 5 runs
# each this is 0.94-1.04 with no sizing field and 3.88-4.13 with one, so 2.5
# separates them with room for fTetWild's run to run variation. An ignored
# field lands at 1.0 and cannot reach it.
MIN_RATIO = 2.5


def _background(fine_below_x: float = 0.0, x_min: float = -0.6) -> pv.UnstructuredGrid:
    """Build a background mesh that asks for fine cells at low x.

    Parameters
    ----------
    fine_below_x : float, default: 0.0
        Points below this x get ``FINE``, the rest get ``COARSE``.
    x_min : float, default: -0.6
        Lower x bound of the background mesh. Raise it to leave part of the
        input uncovered.

    Returns
    -------
    pv.UnstructuredGrid
        All-tetrahedral mesh with a ``target`` point array.

    """
    n = 5
    spacing = (0.6 - x_min) / (n - 1)
    grid = pv.ImageData(
        dimensions=(n, n, n),
        spacing=(spacing, 1.2 / (n - 1), 1.2 / (n - 1)),
        origin=(x_min, -0.6, -0.6),
    )
    bg = grid.to_tetrahedra()
    bg.clear_data()
    bg.point_data["target"] = np.where(bg.points[:, 0] < fine_below_x, FINE, COARSE)
    return bg


def _mean_volume_by_side(grid: pv.UnstructuredGrid) -> tuple[float, float]:
    """Return mean cell volume for x < 0 and x > 0.

    Parameters
    ----------
    grid : pv.UnstructuredGrid
        Tetrahedral mesh to measure.

    Returns
    -------
    tuple[float, float]
        Mean cell volume below and above ``x = 0``.

    """
    sized = grid.compute_cell_sizes(length=False, area=False, volume=True)
    volume = sized.cell_data["Volume"]
    left = grid.cell_centers().points[:, 0] < 0
    assert left.any() and (~left).any(), "both sides must have cells to compare"
    return volume[left].mean(), volume[~left].mean()


def test_sizing_field_refines_the_side_it_asks_for() -> None:
    """A background mesh asking for small cells at low x produces them."""
    mesh = pv.Sphere(theta_resolution=20, phi_resolution=20)
    bg = _background()

    plain_left, plain_right = _mean_volume_by_side(tetrahedralize_pv(mesh, num_opt_iter=5))
    # Without a field the two halves are the same mesh, so this is the
    # baseline the assertion below has to beat.
    assert 0.5 < plain_left / plain_right < 2.0

    sized = tetrahedralize_pv(mesh, num_opt_iter=5, sizing_field=bg, sizing_field_scalars="target")
    left, right = _mean_volume_by_side(sized)
    assert right / left > MIN_RATIO, f"expected the low-x half to be finer, got {left=} {right=}"


def test_sizing_field_leaves_uncovered_regions_alone() -> None:
    """Points outside the background mesh keep the global edge length."""
    mesh = pv.Sphere(theta_resolution=20, phi_resolution=20)
    # Covers x > -0.1 only, and asks for fine cells across all of it, so the
    # uncovered low-x side is the one that must stay coarse.
    bg = _background(fine_below_x=np.inf, x_min=-0.1)
    bg.point_data["target"] = np.full(bg.n_points, FINE)

    left, right = _mean_volume_by_side(
        tetrahedralize_pv(mesh, num_opt_iter=5, sizing_field=bg, sizing_field_scalars="target")
    )
    assert left / right > MIN_RATIO, f"expected the covered half to be finer, got {left=} {right=}"


def test_sizing_field_uses_active_scalars_by_default() -> None:
    """Omitting ``sizing_field_scalars`` picks up the active point array."""
    mesh = pv.Sphere(theta_resolution=20, phi_resolution=20)
    bg = _background()
    bg.set_active_scalars("target")

    left, right = _mean_volume_by_side(tetrahedralize_pv(mesh, num_opt_iter=5, sizing_field=bg))
    assert right / left > MIN_RATIO


def test_sizing_field_low_level_arrays() -> None:
    """``tetrahedralize`` accepts the same field as raw arrays."""
    mesh = pv.Sphere(theta_resolution=20, phi_resolution=20)
    bg = _background()

    points, cells = tetrahedralize(
        mesh.points,
        mesh.faces.reshape(-1, 4)[:, 1:],
        num_opt_iter=5,
        vtk_ordering=True,
        bg_vertices=bg.points,
        bg_tets=bg.cells_dict[pv.CellType.TETRA],
        bg_values=bg.point_data["target"],
    )
    grid = pv.UnstructuredGrid({pv.CellType.TETRA: cells}, points)
    left, right = _mean_volume_by_side(grid)
    assert right / left > MIN_RATIO


def test_sizing_field_requires_optimize() -> None:
    """fTetWild only applies the field while optimizing, so refuse otherwise."""
    mesh = pv.Sphere(theta_resolution=10, phi_resolution=10)
    with pytest.raises(ValueError, match="requires `optimize=True`"):
        tetrahedralize_pv(
            mesh, optimize=False, sizing_field=_background(), sizing_field_scalars="target"
        )


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"bg_vertices": np.zeros((4, 3))}, ValueError, "must all be given"),
        (
            {
                "bg_vertices": np.zeros((4, 3)),
                "bg_tets": np.zeros((1, 4), np.int32),
                "bg_values": np.zeros(3),
            },
            ValueError,
            "one value per point",
        ),
        (
            {
                "bg_vertices": np.zeros((4, 3)),
                "bg_tets": np.full((1, 4), 9, np.int32),
                "bg_values": np.ones(4),
            },
            IndexError,
            "outside `bg_vertices`",
        ),
        (
            {
                "bg_vertices": np.zeros((4, 3)),
                "bg_tets": np.zeros((1, 4), np.int32),
                "bg_values": np.zeros(4),
            },
            ValueError,
            "finite and positive",
        ),
    ],
)
def test_sizing_field_array_validation(kwargs: dict, error: type, match: str) -> None:
    """Malformed low-level arrays are rejected before reaching the extension."""
    mesh = pv.Sphere(theta_resolution=10, phi_resolution=10)
    with pytest.raises(error, match=match):
        tetrahedralize(mesh.points, mesh.faces.reshape(-1, 4)[:, 1:], **kwargs)


def test_sizing_field_rejects_non_tetrahedral_grid() -> None:
    """The background mesh has to be tetrahedra, which is what fTetWild reads."""
    mesh = pv.Sphere(theta_resolution=10, phi_resolution=10)
    hexes = pv.ImageData(dimensions=(3, 3, 3)).cast_to_unstructured_grid()
    hexes.point_data["target"] = np.full(hexes.n_points, FINE)
    with pytest.raises(ValueError, match="only tetrahedra"):
        tetrahedralize_pv(mesh, sizing_field=hexes, sizing_field_scalars="target")


def test_sizing_field_rejects_wrong_type() -> None:
    """A surface mesh is a common mistake and gets a clear message."""
    mesh = pv.Sphere(theta_resolution=10, phi_resolution=10)
    with pytest.raises(TypeError, match="must be a pyvista.UnstructuredGrid"):
        tetrahedralize_pv(mesh, sizing_field=pv.Sphere())


def test_sizing_field_missing_scalars() -> None:
    """A named array that is not there is an error, not a silent default."""
    mesh = pv.Sphere(theta_resolution=10, phi_resolution=10)
    with pytest.raises(KeyError, match="no point array named 'nope'"):
        tetrahedralize_pv(mesh, sizing_field=_background(), sizing_field_scalars="nope")


def test_sizing_field_no_active_scalars() -> None:
    """With nothing active and no name given, say which argument to pass."""
    mesh = pv.Sphere(theta_resolution=10, phi_resolution=10)
    bg = _background()
    bg.set_active_scalars(None)
    with pytest.raises(ValueError, match="no active point scalars"):
        tetrahedralize_pv(mesh, sizing_field=bg)
