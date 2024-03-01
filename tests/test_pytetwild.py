import numpy as np
import pyvista as pv
import pytest
from your_module import tetrahedralize_pv, tetrahedralize  # Replace 'your_module' with the name of your Python file

# Parameterized test for tetrahedralize_pv function
@pytest.mark.parametrize("mesh_generator", [pv.Icosphere, pv.examples.download_bunny_coarse])
def test_tetrahedralize_pv(mesh_generator):
    mesh = mesh_generator()
    result = tetrahedralize_pv(mesh)
    assert isinstance(result, pv.UnstructuredGrid), "The result should be a PyVista UnstructuredGrid"
    assert result.n_cells > 0, "The resulting mesh should have more than 0 cells"
    assert result.n_points > 0, "The resulting mesh should have more than 0 points"

# Parameterized test for tetrahedralize function
@pytest.mark.parametrize("mesh_generator", [pv.Icosphere, pv.examples.download_bunny_coarse])
def test_tetrahedralize(mesh_generator):
    mesh = mesh_generator()
    vertices = np.array(mesh.points, dtype=np.float64)
    faces = np.array(mesh.faces.reshape((-1, 4))[:, 1:4], dtype=np.int32)
    
    vertices_result, tetrahedra_result = tetrahedralize(vertices, faces)
    assert isinstance(vertices_result, np.ndarray), "The vertices result should be a numpy array"
    assert isinstance(tetrahedra_result, np.ndarray), "The tetrahedra result should be a numpy array"
    assert len(vertices_result) > 0, "There should be more than 0 vertices in the result"
    assert len(tetrahedra_result) > 0, "There should be more than 0 tetrahedra in the result"
