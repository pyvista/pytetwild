"""Wrapper for fTetWild."""

import warnings
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from pytetwild import PyfTetWildWrapper

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


def _unpack_sizing_field(
    sizing_field: "pv.UnstructuredGrid | None",
    scalars: str | None,
) -> tuple[NDArray[np.float64] | None, NDArray[np.int32] | None, NDArray[np.float64] | None]:
    """Split a PyVista sizing field into points, tetrahedra, and values.

    Parameters
    ----------
    sizing_field : pv.UnstructuredGrid | None
        All-tetrahedral background mesh, or ``None`` for no sizing field.
    scalars : str | None
        Name of the point array holding the target edge length. ``None``
        uses the active point scalars.

    Returns
    -------
    tuple
        ``(points, tets, values)``, each ``None`` when ``sizing_field`` is
        ``None``. Validation of the arrays themselves is left to
        :func:`_check_sizing_field`.

    """
    if sizing_field is None:
        return None, None, None

    import pyvista.core as pv

    if not isinstance(sizing_field, pv.UnstructuredGrid):
        raise TypeError(
            f"`sizing_field` must be a pyvista.UnstructuredGrid, got {type(sizing_field)}"
        )

    cells_by_type = sizing_field.cells_dict
    if list(cells_by_type) != [VTK_TETRA]:
        raise ValueError(
            "`sizing_field` must contain only tetrahedra. Call `.clean()` on a"
            " tetrahedralized mesh or build one with e.g."
            " `pyvista.ImageData(...).to_tetrahedra()`."
        )

    if scalars is None:
        scalars = sizing_field.point_data.active_scalars_name
        if not scalars:
            raise ValueError(
                "`sizing_field` has no active point scalars. Pass `sizing_field_scalars`"
                " with the name of the target edge length array."
            )
    if scalars not in sizing_field.point_data:
        raise KeyError(f"`sizing_field` has no point array named {scalars!r}")

    return sizing_field.points, cells_by_type[VTK_TETRA], sizing_field.point_data[scalars]


def _check_sizing_field(
    bg_vertices: NDArray[np.float64] | None,
    bg_tets: NDArray[np.int32] | None,
    bg_values: NDArray[np.float64] | None,
    optimize: bool,
) -> tuple[NDArray[np.float64], NDArray[np.int32], NDArray[np.float64]]:
    """Validate a sizing field and return it as the arrays the extension takes.

    Parameters
    ----------
    bg_vertices : np.ndarray[np.float64] | None
        Background mesh points, shape ``(n, 3)``.
    bg_tets : np.ndarray[np.int32] | None
        Background mesh tetrahedra, shape ``(m, 4)``.
    bg_values : np.ndarray[np.float64] | None
        Target edge length at each background point, shape ``(n,)``.
    optimize : bool
        Whether optimization is enabled. fTetWild applies the sizing field
        from its optimization stage, so a field is silently ignored without
        it and that is worth an error rather than a surprise.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        Vertices, tetrahedra, and values. All three are empty when no
        sizing field was given.

    """
    given = [arr is not None for arr in (bg_vertices, bg_tets, bg_values)]
    if not any(given):
        return (
            np.empty((0, 3), dtype=np.float64),
            np.empty((0, 4), dtype=np.int32),
            np.empty(0, dtype=np.float64),
        )
    if not all(given):
        raise ValueError(
            "`bg_vertices`, `bg_tets`, and `bg_values` must all be given to use a sizing field"
        )
    if not optimize:
        raise ValueError(
            "A sizing field requires `optimize=True`. fTetWild applies it during "
            "optimization, so it would have no effect."
        )

    vertices = np.ascontiguousarray(bg_vertices, dtype=np.float64)
    tets = np.ascontiguousarray(bg_tets, dtype=np.int32)
    values = np.ascontiguousarray(bg_values, dtype=np.float64).ravel()

    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError(f"`bg_vertices` must have shape (n, 3), got {vertices.shape}")
    if tets.ndim != 2 or tets.shape[1] != 4:
        raise ValueError(f"`bg_tets` must have shape (m, 4), got {tets.shape}")
    if len(vertices) == 0 or len(tets) == 0:
        raise ValueError("Sizing field must have at least one point and one tetrahedron")
    if len(values) != len(vertices):
        raise ValueError(
            f"`bg_values` must have one value per point, got {len(values)} "
            f"for {len(vertices)} points"
        )
    if tets.min() < 0 or tets.max() >= len(vertices):
        raise IndexError(
            f"`bg_tets` indexes points outside `bg_vertices` (0 to {len(vertices) - 1})"
        )
    if not np.isfinite(vertices).all():
        raise ValueError("`bg_vertices` must be finite")
    if not np.isfinite(values).all() or values.min() <= 0:
        raise ValueError("`bg_values` are target edge lengths and must be finite and positive")

    return vertices, tets, values


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

        # PyVista chooses the VTK build -- stock vtkmodules or the cvista fork
        # -- and arrays from the other one are a foreign type to the grid they
        # are handed to, so take these names from whichever it selected.
        from pyvista._vtk import numpy_to_vtk, vtkCellArray, vtkTypeInt32Array
    except ImportError as exc:
        raise ModuleNotFoundError(
            "Install PyVista to use this feature with:\n\npip install pytetwild[all]"
        ) from exc

    if cells.ndim != 2:
        raise ValueError("cells must be 2D")

    grid = pv.UnstructuredGrid()
    grid.points = points

    n_cells, n_nodes_per_cell = cells.shape
    vtk_dtype = vtkTypeInt32Array().GetDataType()
    conn_vtk = numpy_to_vtk(cells.ravel(), deep=False, array_type=vtk_dtype)

    # Every cell here is the same width, so VTK 9.6.2 can hold that one width
    # instead of an offset per cell and synthesize the offsets on demand. That
    # is an (n_cells + 1) int32 array we no longer build or keep.
    cell_array = vtkCellArray()
    if not cell_array.SetData(n_nodes_per_cell, conn_vtk):
        offsets = np.arange(0, n_cells * n_nodes_per_cell + 1, n_nodes_per_cell, dtype=np.int32)
        cell_array.SetData(numpy_to_vtk(offsets, deep=False, array_type=vtk_dtype), conn_vtk)

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
    quiet: bool = True,
    disable_filtering: bool = False,
    sizing_field: "pv.UnstructuredGrid | None" = None,
    sizing_field_scalars: str | None = None,
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
    disable_filtering : bool, default: False
        Disable the filtering of the resulting mesh, and thus keep
        fTetWilds background mesh.
    sizing_field : pv.UnstructuredGrid, optional
        All-tetrahedral mesh supplying a spatially varying target edge
        length, so some regions come out finer than others. The value at any
        point is interpolated over the containing tetrahedron; anywhere this
        mesh does not cover falls back to the global ideal edge length.
        Requires ``optimize=True``. This is fTetWild's ``--bg-mesh``, and is
        unrelated to ``disable_filtering``, which keeps the internal
        background mesh fTetWild builds for itself.
    sizing_field_scalars : str, optional
        Name of the point array on ``sizing_field`` holding the target edge
        length, in the same units as the input. Defaults to its active point
        scalars.

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
    bg_v, bg_t, bg_val = _check_sizing_field(
        *_unpack_sizing_field(sizing_field, sizing_field_scalars), optimize=optimize
    )

    if edge_length_abs is None:
        edge_length_abs = 0.0

    vertices = mesh.points.astype(np.float64, copy=False)
    faces = mesh._connectivity_array.reshape(-1, 3)
    if faces.dtype == np.int32:
        # we can cheat here and just treat it as unsigned 32 (assuming no negative indices)
        faces_unsigned = faces.view(np.uint32)
    else:
        faces_unsigned = faces.astype(np.uint32, copy=False)

    skip_simplify = not simplify
    vtk_ordering = True
    tmesh_v, tmesh_c = PyfTetWildWrapper.tetrahedralize_mesh(
        vertices,
        faces_unsigned,
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
        disable_filtering,
        bg_v,
        bg_t,
        bg_val,
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
    quiet: bool = True,
    vtk_ordering: bool = False,
    disable_filtering: bool = False,
    bg_vertices: NDArray[np.float64] | None = None,
    bg_tets: NDArray[np.int32] | None = None,
    bg_values: NDArray[np.float64] | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.int32]]:
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
    vtk_ordering : bool, default: False
        Reorder the tetrahedral cell indices to match VTK's ordering.
    disable_filtering : bool, default: False
        Disable the filtering of the resulting mesh, and thus keep
        fTetWilds background mesh.
    bg_vertices : np.ndarray[np.float64], optional
        Points of a background tetrahedral mesh supplying a spatially
        varying target edge length, shape ``(n, 3)``. Give all three
        ``bg_`` arguments together, and see ``bg_values``.
    bg_tets : np.ndarray[np.int32], optional
        Tetrahedra of the background mesh, shape ``(m, 4)``, indexing
        ``bg_vertices``.
    bg_values : np.ndarray[np.float64], optional
        Target edge length at each background point, shape ``(n,)``, in the
        same units as the input. The value at any point is interpolated over
        the containing background tetrahedron; anywhere the background mesh
        does not cover falls back to the global ideal edge length. Requires
        ``optimize=True``.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing the vertices and tetrahedra of the tetrahedral mesh.
    """
    _check_edge_length(edge_length_fac)
    bg_v, bg_t, bg_val = _check_sizing_field(bg_vertices, bg_tets, bg_values, optimize=optimize)
    if not isinstance(vertices, np.ndarray):
        raise TypeError("`vertices` must be a numpy array")
    if not isinstance(faces, np.ndarray):
        raise TypeError("`faces` must be a numpy array")
    vertices = vertices.astype(np.float64, copy=False)
    if faces.dtype == np.int32:
        # we can cheat here and just treat it as unsigned 32 (assuming no negative indices)
        faces_unsigned = faces.view(np.uint32)
    else:
        faces_unsigned = faces.astype(np.uint32, copy=False)

    if edge_length_abs is None:
        edge_length_abs = 0.0

    skip_simplify = not simplify
    return PyfTetWildWrapper.tetrahedralize_mesh(
        vertices,
        faces_unsigned,
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
        disable_filtering,
        bg_v,
        bg_t,
        bg_val,
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
