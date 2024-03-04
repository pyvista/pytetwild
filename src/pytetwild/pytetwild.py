"""Wrapper for fTetWild."""
import warnings
import numpy as np
from pytetwild import PyfTetWildWrapper
from typing import Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    import pyvista as pv


def tetrahedralize_pv(mesh: "pv.PolyData") -> "pv.UnstructuredGrid":
    """
    Convert a PyVista surface mesh to a PyVista unstructured grid.

    Parameters
    ----------
    mesh : pv.PolyData
        The input surface mesh.

    Returns
    -------
    pv.UnstructuredGrid
        The converted unstructured grid.
    """
    try:
        import pyvista as pv
    except:
        raise ModuleNotFoundError(
            "Install PyVista to use this feature with:\n\n" "pip install pytetwild[all]"
        )

    if not mesh.is_all_triangles:
        warnings.warn(
            "Input mesh is not all triangles. Either call `.triangulate()`"
            " beforehand to suppress this warning or use an all triangle mesh."
        )
        mesh = mesh.triangulate()

    vertices = np.array(mesh.points, dtype=np.float64)
    faces = np.array(mesh.faces.reshape((-1, 4))[:, 1:4], dtype=np.int32)

    (
        tetrahedral_mesh_vertices,
        tetrahedral_mesh_tetrahedra,
    ) = PyfTetWildWrapper.tetrahedralize_mesh(vertices, faces)

    cells = np.hstack(
        [
            np.full((tetrahedral_mesh_tetrahedra.shape[0], 1), 4, dtype=np.int32),
            tetrahedral_mesh_tetrahedra,
        ]
    )
    cell_types = np.full(tetrahedral_mesh_tetrahedra.shape[0], 10, dtype=np.uint8)

    return pv.UnstructuredGrid(cells, cell_types, tetrahedral_mesh_vertices)


def tetrahedralize(vertices: np.ndarray, faces: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert mesh vertices and faces to a tetrahedral mesh.

    Parameters
    ----------
    vertices : np.ndarray[double]
        The vertices of the mesh.
    faces : np.ndarray
        The faces of the mesh.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing the vertices and tetrahedra of the tetrahedral mesh.
    """
    return PyfTetWildWrapper.tetrahedralize_mesh(vertices, faces)
