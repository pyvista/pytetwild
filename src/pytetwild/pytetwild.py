import numpy as np
import pyvista as pv
from pytetwild import PyfTetWildWrapper
from typing import Tuple

def tetrahedralize_pv(mesh: pv.PolyData) -> pv.UnstructuredGrid:
    """
    Converts a PyVista surface mesh to a PyVista unstructured grid.

    Parameters
    ----------
    mesh : pv.PolyData
        The input surface mesh.

    Returns
    -------
    pv.UnstructuredGrid
        The converted unstructured grid.
    """
    vertices = np.array(mesh.points, dtype=np.float64)
    faces = np.array(mesh.faces.reshape((-1, 4))[:, 1:4], dtype=np.int32)
    
    tetrahedral_mesh_vertices, tetrahedral_mesh_tetrahedra = PyfTetWildWrapper.tetrahedralize_mesh(vertices, faces)
    
    cells = np.hstack([np.full((tetrahedral_mesh_tetrahedra.shape[0], 1), 4, dtype=np.int32), tetrahedral_mesh_tetrahedra])
    cell_types = np.full(tetrahedral_mesh_tetrahedra.shape[0], 10, dtype=np.uint8)
    
    return pv.UnstructuredGrid(cells, cell_types, tetrahedral_mesh_vertices)

def tetrahedralize(vertices: np.ndarray, faces: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Converts mesh vertices and faces to a tetrahedral mesh using PyfTetWildWrapper.

    Parameters
    ----------
    vertices : np.ndarray
        The vertices of the mesh.
    faces : np.ndarray
        The faces of the mesh.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing the vertices and tetrahedra of the tetrahedral mesh.
    """
    return PyfTetWildWrapper.tetrahedralize_mesh(vertices, faces)
