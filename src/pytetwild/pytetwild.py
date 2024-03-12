"""Wrapper for fTetWild."""
import warnings
import numpy as np
from pytetwild import PyfTetWildWrapper  # type: ignore
from typing import Tuple, TYPE_CHECKING
import numpy.typing as npt

if TYPE_CHECKING:
    import pyvista as pv


def _check_edge_length(edge_length_fac: float) -> None:
    """Check edge length.

    Parameters
    ----------
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagional. The
        default ideal edge length is bb/20 (bounding box divided by 20).


    """
    if edge_length_fac > 1.0 or edge_length_fac < 0.0 + 1e-16:
        raise ValueError("Edge length factor must be between 1e-16 and 1.0")


def tetrahedralize_pv(
    mesh: "pv.PolyData", edge_length_fac: float = 0.05, optimize: bool = True
) -> "pv.UnstructuredGrid":
    """
    Convert a PyVista surface mesh to a PyVista unstructured grid.

    Parameters
    ----------
    mesh : pv.PolyData
        The input surface mesh.
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagional. The
        default ideal edge length is bb/20 (bounding box divided by 20).
    optimize : bool
        Improve the minimum scaled Jacobean for each cell. This leads to higher
        cell quality at the expense of computation time.

    Returns
    -------
    pv.UnstructuredGrid
        The converted unstructured grid containing only tetrahedra.

    Examples
    --------
    >>> import pyvista as pv
    >>> import pytetwild
    >>> surface_mesh = pv.Sphere()
    >>> tetrahedral_mesh = pytetwild.tetrahedralize_pv(surface_mesh)
    >>> tetrahedral_mesh
    UnstructuredGrid (0x7ff568593d00)
      N Cells:    7247
      N Points:   1658
      X Bounds:   -4.993e-01, 4.993e-01
      Y Bounds:   -4.965e-01, 4.965e-01
      Z Bounds:   -5.000e-01, 4.999e-01
      N Arrays:   0

    """
    try:
        import pyvista as pv
    except:
        raise ModuleNotFoundError(
            "Install PyVista to use this feature with:\n\n" "pip install pytetwild[all]"
        )

    if not isinstance(mesh, pv.PolyData):
        raise TypeError(f"`mesh` must be a pyvista.PolyData, got {type(mesh)}")

    if not mesh.is_all_triangles:
        warnings.warn(
            "Input mesh is not all triangles. Either call `.triangulate()`"
            " beforehand to suppress this warning or use an all triangle mesh."
        )
        mesh = mesh.triangulate()

    _check_edge_length(edge_length_fac)
    vertices = mesh.points.astype(np.float64, copy=False)
    faces = mesh.faces.reshape((-1, 4))[:, 1:4].astype(np.float32, copy=False)

    (
        tetrahedral_mesh_vertices,
        tetrahedral_mesh_tetrahedra,
    ) = PyfTetWildWrapper.tetrahedralize_mesh(vertices, faces, optimize, edge_length_fac)

    cells = np.hstack(
        [
            np.full((tetrahedral_mesh_tetrahedra.shape[0], 1), 4, dtype=np.int32),
            tetrahedral_mesh_tetrahedra,
        ]
    )
    cell_types = np.full(tetrahedral_mesh_tetrahedra.shape[0], 10, dtype=np.uint8)

    return pv.UnstructuredGrid(cells, cell_types, tetrahedral_mesh_vertices)


def tetrahedralize(
    vertices: npt.NDArray[np.float64],
    faces: npt.NDArray[np.int32],
    optimize: bool = True,
    edge_length_fac: float = 0.05,
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.int32]]:
    """
    Convert mesh vertices and faces to a tetrahedral mesh.

    Parameters
    ----------
    vertices : np.ndarray[double]
        The vertices of the mesh.
    faces : np.ndarray
        The faces of the mesh.
    optimize : bool
        Improve the minimum scaled Jacobean for each cell. This leads to higher
        cell quality at the expense of computation time.
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagional. The
        default ideal edge length is bb/20 (bounding box divided by 20).

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing the vertices and tetrahedra of the tetrahedral mesh.
    """
    _check_edge_length(edge_length_fac)
    if not isinstance(vertices, np.ndarray):
        raise TypeError("`vertices` must be a numpy array")
    if not isinstance(faces, np.ndarray):
        raise TypeError("`faces` must be a numpy array")
    vertices = vertices.astype(np.float64, copy=False)
    faces = faces.astype(np.int32, copy=False)
    return PyfTetWildWrapper.tetrahedralize_mesh(vertices, faces, optimize, edge_length_fac)
