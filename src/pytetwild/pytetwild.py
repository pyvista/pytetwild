"""Wrapper for fTetWild."""

import warnings
import numpy as np
from pytetwild import PyfTetWildWrapper
from typing import TYPE_CHECKING
from numpy.typing import NDArray

if TYPE_CHECKING:
    import pyvista as pv

VTK_UNSIGNED_CHAR = 3
VTK_TETRA = 10
VTK_QUADRATIC_TETRA = 24


def _check_edge_length(edge_length_fac: float) -> None:
    """Check edge length.

    Parameters
    ----------
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagonal. The
        default ideal edge length is bb/20 (bounding box divided by 20).

    """
    if edge_length_fac > 1.0 or edge_length_fac < 0.0 + 1e-16:
        raise ValueError("Edge length factor must be between 1e-16 and 1.0")


def _ugrid_from_regular_cells(
    points: NDArray[np.float32] | NDArray[np.float64],
    cells: NDArray[np.int32],
) -> "pv.UnstructuredGrid":
    """
    Create an UnstructuredGrid from points and fixed-size cells.

    Parameters
    ----------
    points : np.ndarray[np.float64]
        Point coordinates.
    cells : np.ndarray[np.int32]
        Cell connectivity without padding. Cell type inferred from shape.

    Returns
    -------
    pyvista.UnstructuredGrid
        Unstructured grid.

    """
    try:
        import pyvista.core as pv
    except:
        raise ModuleNotFoundError(
            "Install PyVista to use this feature with:\n\npip install pytetwild[all]"
        )

    from vtkmodules.vtkCommonCore import vtkTypeInt32Array
    from vtkmodules.util.numpy_support import numpy_to_vtk
    from vtkmodules.vtkCommonDataModel import vtkCellArray

    if cells.ndim != 2:
        raise ValueError("cells must be 2D")

    grid = pv.UnstructuredGrid()
    grid.points = points

    n_cells, n_nodes_per_cell = cells.shape
    vtk_dtype = vtkTypeInt32Array().GetDataType()
    offsets = np.arange(0, n_cells * n_nodes_per_cell + 1, n_nodes_per_cell, dtype=np.int32)
    offsets_vtk = numpy_to_vtk(offsets, deep=False, array_type=vtk_dtype)
    conn_vtk = numpy_to_vtk(cells.ravel(), deep=False, array_type=vtk_dtype)

    cell_array = vtkCellArray()
    cell_array.SetData(offsets_vtk, conn_vtk)

    if n_nodes_per_cell == 4:
        cell_type = VTK_TETRA
    elif n_nodes_per_cell == 10:
        cell_type = VTK_QUADRATIC_TETRA
    else:
        raise ValueError(f"Unsupported number of nodes per cells {n_nodes_per_cell}")

    celltypes = numpy_to_vtk(
        np.full(n_cells, cell_type, dtype=np.uint8), deep=False, array_type=VTK_UNSIGNED_CHAR
    )
    grid.SetCells(celltypes, cell_array)
    return grid


def tetrahedralize_pv(
    mesh: "pv.PolyData",
    edge_length_fac: float = 0.05,
    edge_length_abs: float | None = None,
    optimize: bool = True,
    simplify: bool = True,
    epsilon: float = 1e-3,
    stop_energy: float = 10.0,
    coarsen: bool = False,
    num_threads: int = 0,
    num_opt_iter: int = 80,
    loglevel: int = 3,
    quiet: bool = False,
) -> "pv.UnstructuredGrid":
    """
    Convert a PyVista surface mesh to a PyVista unstructured grid.

    Parameters
    ----------
    mesh : pv.PolyData
        The input surface mesh. Should be composed of all triangles.
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagonal. The
        default ideal edge length is ``bb/20`` (bounding box divided by
        20). Ignored when ``edge_length_abs`` is input.
    edge_length_abs : float, optional
        Absolute ideal edge length. When input ``edge_length_fac`` is ignored.
    optimize : bool, default: True
        Improve the minimum scaled Jacobean for each cell. This leads to higher
        cell quality at the expense of computation time. Optimization level is
        dependent on ``stop_energy`` and ``num_opt_iter``.
    simplify : bool, default: True
        Simplfiy the input mesh surface before tetrahedralization.
    epsilon : float, default 1e-3
        Envelop size, specifying the maximum distance of the output surface
        from the input surface, relative to the bounding box size.
    stop_energy : float, default: 10.0
        The mesh optimization stops when the conformal AMIPS energy reaches
        ``stop_energy``.
    coarsen : bool, default: False
        Coarsen the output as much as possible, while maintaining the mesh
        quality.
    num_threads : int, default: 0
        Set number of threads used. 0 (default) uses all available cores.
    num_opt_iter : int, default: 80
        Maximum number of optimization iterations if ``optimize=True``.
    loglevel : int, default: 6
        Set log level (0 = most verbose, 6 = minimal output).
    quiet : bool, default: False
        Disable all output. Overrides ``loglevel``.

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
        import pyvista.core as pv
    except:
        raise ModuleNotFoundError(
            "Install PyVista to use this feature with:\n\npip install pytetwild[all]"
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

    if edge_length_abs is None:
        edge_length_abs = 0.0

    vertices = mesh.points
    faces = mesh._connectivity_array.reshape(-1, 3)
    skip_simplify = not simplify
    vtk_ordering = True
    tmesh_v, tmesh_c = PyfTetWildWrapper.tetrahedralize_mesh(
        vertices,
        faces,
        optimize,
        skip_simplify,
        edge_length_fac,
        edge_length_abs,
        epsilon,
        stop_energy,
        coarsen,
        num_threads,
        num_opt_iter,
        loglevel,
        quiet,
        vtk_ordering,
    )
    return _ugrid_from_regular_cells(tmesh_v, tmesh_c)


def tetrahedralize(
    vertices: NDArray[np.float32] | NDArray[np.float64],
    faces: NDArray[np.int32],
    edge_length_fac: float = 0.05,
    edge_length_abs: float | None = None,
    optimize: bool = True,
    simplify: bool = True,
    epsilon: float = 1e-3,
    stop_energy: float = 10.0,
    coarsen: bool = False,
    num_threads: int = 0,
    num_opt_iter: int = 80,
    loglevel: int = 3,
    quiet: bool = False,
) -> tuple[NDArray[np.float32], NDArray[np.int32]]:
    """
    Convert mesh vertices and faces to a tetrahedral mesh.

    Parameters
    ----------
    vertices : np.ndarray[np.float32] | np.ndarray[np.float64]
        The vertices of the mesh.
    faces : np.ndarray[np.int32]
        The faces of the mesh.
    edge_length_fac : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagonal. The
        default ideal edge length is bb/20 (bounding box divided by
        20). Ignored when ``edge_length_abs`` is input.
    edge_length_abs : float, optional
        Absolute ideal edge length. When input ``edge_length_fac`` is ignored.
    optimize : bool
        Improve the minimum scaled Jacobean for each cell. This leads to higher
        cell quality at the expense of computation time.
    simplify : bool, default: True
        Simplfiy the input mesh surface before tetrahedralization.
    epsilon : float, default 1e-3
        Envelop size, specifying the maximum distance of the output surface
        from the input surface, relative to the bounding box size.
    stop_energy : float, default: 10.0
        The mesh optimization stops when the conformal AMIPS energy reaches
        ``stop_energy``.
    coarsen : bool, default: False
        Coarsen the output as much as possible, while maintaining the mesh
        quality.
    num_threads : int, default: 0
        Set number of threads used. 0 (default) uses all available cores.
    num_opt_iter : int, default: 80
        Maximum number of optimization iterations if ``optimize=True``.
    loglevel : int, default: 6
        Set log level (0 = most verbose, 6 = minimal output).
    quiet : bool, default: False
        Disable all output. Overrides ``loglevel``.

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

    if edge_length_abs is None:
        edge_length_abs = 0.0

    skip_simplify = not simplify
    vtk_ordering = False
    return PyfTetWildWrapper.tetrahedralize_mesh(
        vertices,
        faces,
        optimize,
        skip_simplify,
        edge_length_fac,
        edge_length_abs,
        epsilon,
        stop_energy,
        coarsen,
        num_threads,
        num_opt_iter,
        loglevel,
        quiet,
        vtk_ordering,
    )


def tetrahedralize_csg(
    csg_file: str,
    epsilon: float = 1e-3,
    edge_length_r: float = 0.05,
    stop_energy: float = 10.0,
    coarsen: bool = True,
    num_threads: int = 0,
    loglevel: int = 3,
) -> "pv.UnstructuredGrid":
    """
    Generate a tetrahedral mesh based on a the CSG tree specified in the csf_file.

    Parameters
    ----------
    csg_file : str
        Path to the input json file.
    epsilon : float, default 1e-3
        Envelop size, specifying the maximum distance of the output surface
        from the input surface, relative to the bounding box size.
    edge_length_r : float, default: 0.05
        Tetrahedral edge length as a function of bounding box diagonal. The
        default ideal edge length is bb/20 (bounding box divided by 20).
    stop_energy : float, default: 10.0
        The mesh optimization stops when the  conformal AMIPS energy reaches 'stop_energy'.
    coarsen : bool, default: true
        Coarsen the output as much as possible, while maintaining the mesh quality.
    num_threads : int, default: 0
        Set number of threads used (0 means all available cores).
    loglevel : int, default: 6
        Set log level (0 = most verbose, 6 = minimal output).

    Returns
    -------
    pv.UnstructuredGrid
        The converted unstructured grid containing only tetrahedra, with a cell
        attribute 'marker' indicating which of the input surfaces the cell
        belongs to.
    """
    try:
        pass
    except:
        raise ModuleNotFoundError(
            "Install PyVista to use this feature with:\n\npip install pytetwild[all]"
        )
    vtk_ordering = True
    tmesh_v, tmesh_c, tmesh_marker = PyfTetWildWrapper.tetrahedralize_csg(
        csg_file,
        epsilon,
        edge_length_r,
        stop_energy,
        coarsen,
        num_threads,
        loglevel,
        vtk_ordering,
    )

    ugrid = _ugrid_from_regular_cells(tmesh_v, tmesh_c).clean()
    ugrid["marker"] = tmesh_marker
    return ugrid
