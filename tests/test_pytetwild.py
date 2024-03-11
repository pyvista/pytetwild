import os
import numpy as np
import pyvista as pv
import vtk
from scipy.spatial import KDTree
import pytest
from pytetwild import (
    tetrahedralize_pv,
    tetrahedralize,
)

THIS_PATH = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def default_test_data():
    data = {
        "input": pv.read(os.path.join(THIS_PATH, "test_data/test_surf.ply")),
        "output": pv.read(os.path.join(THIS_PATH, "test_data/test_tets.msh")),
    }
    return data


# Parameterized test for tetrahedralize_pv function
@pytest.mark.parametrize("mesh_generator", [pv.Icosphere, pv.examples.download_bunny_coarse])
def test_tetrahedralize_pv(mesh_generator):
    mesh = mesh_generator()
    result = tetrahedralize_pv(mesh, edge_length_fac=0.5)
    assert isinstance(
        result, pv.UnstructuredGrid
    ), "The result should be a PyVista UnstructuredGrid"
    assert result.n_cells > 0, "The resulting mesh should have more than 0 cells"
    assert result.n_points > 0, "The resulting mesh should have more than 0 points"


def test_tetrahedralize_edge_length():
    mesh = pv.Cube().triangulate()
    result = tetrahedralize_pv(mesh)
    result_very_coarse = tetrahedralize_pv(mesh, edge_length_fac=1.0)
    assert result_very_coarse.n_cells < result.n_cells
    with pytest.raises(ValueError):
        tetrahedralize_pv(mesh, edge_length_fac=0.0)


def test_tetrahedralize_pv_opt():
    mesh = pv.Sphere(phi_resolution=10, theta_resolution=10)
    grid = tetrahedralize_pv(mesh, optimize=True)
    qual_mean = np.mean(-grid.compute_cell_quality()["CellQuality"])

    grid_no_opt = tetrahedralize_pv(mesh, optimize=False)
    qual_mean_no_opt = np.mean(-grid_no_opt.compute_cell_quality()["CellQuality"])
    assert qual_mean > qual_mean_no_opt


# Parameterized test for tetrahedralize function
@pytest.mark.parametrize("mesh_generator", [pv.Icosphere, pv.examples.download_bunny_coarse])
def test_tetrahedralize(mesh_generator):
    mesh = mesh_generator()
    vertices = mesh.points
    faces = mesh.faces.reshape((-1, 4))[:, 1:4]

    vertices_result, tetrahedra_result = tetrahedralize(vertices, faces, edge_length_fac=0.5)
    assert isinstance(vertices_result, np.ndarray), "The vertices result should be a numpy array"
    assert isinstance(
        tetrahedra_result, np.ndarray
    ), "The tetrahedra result should be a numpy array"
    assert len(vertices_result) > 0, "There should be more than 0 vertices in the result"
    assert len(tetrahedra_result) > 0, "There should be more than 0 tetrahedra in the result"


def _sample_points_vtk(mesh_pv, dist_btw_pts=0.01):
    point_sampler = vtk.vtkPolyDataPointSampler()
    point_sampler.SetInputData(mesh_pv)
    point_sampler.SetDistance(dist_btw_pts)
    point_sampler.SetPointGenerationMode(point_sampler.REGULAR_GENERATION)
    point_sampler.Update()
    points_sampled = pv.PolyData(point_sampler.GetOutput()).points
    return points_sampled


def _symmetric_surf_dist(pts0, pts1):
    d_kdtree0, _ = KDTree(pts0).query(pts1)
    d_kdtree1, _ = KDTree(pts1).query(pts0)
    return (np.mean(d_kdtree0) + np.mean(d_kdtree1)) / 2


@pytest.mark.parametrize(
    "mesh_generator", [pv.Icosphere]
)  # pv.examples.download_bunny_coarse is not closed, so select_enclosed_points fails
def test_output_points_enclosed(mesh_generator):
    input_pv = mesh_generator()
    py_output_pv = tetrahedralize_pv(input_pv, edge_length_fac=0.1)
    additional_input_scaling = 0.01
    enclosed_pv = py_output_pv.select_enclosed_points(input_pv.scale(1 + additional_input_scaling))
    enclosed_ratio = enclosed_pv.point_data["SelectedPoints"].sum() / input_pv.points.shape[0]
    assert (
        enclosed_ratio > 0.99
    ), "all output vertices should be within some threshold of the input surf"


def test_default_output_surf_dist(default_test_data):
    input_pv = default_test_data["input"]
    output_pv = default_test_data["output"]
    py_output_pv = tetrahedralize_pv(input_pv)
    pts0 = _sample_points_vtk(py_output_pv.extract_surface())
    pts1 = _sample_points_vtk(output_pv.extract_surface())
    surf_dist = _symmetric_surf_dist(pts0, pts1)
    assert surf_dist < 1e-2, "surfs of outputs from c++/py should be similar"
