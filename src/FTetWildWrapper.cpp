#include <floattetwild/AABBWrapper.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/Types.hpp>

#include <memory>
#include <oneapi/tbb/global_control.h>
#include <thread>

#include <Eigen/Dense>
#include <iostream>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <geogram/api/defs.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>

#include "array_support.h"
#include "igl/default_num_threads.h"

namespace nb = nanobind;
using namespace nb::literals;

using floatTetWild::AABBWrapper;
using floatTetWild::boolean_operation;
using floatTetWild::CSGTreeParser;
using floatTetWild::FloatTetDelaunay;
using floatTetWild::json;
using floatTetWild::Vector3;
using floatTetWild::Vector3i;

struct CoutRedirect {
    std::streambuf *old = nullptr;
    explicit CoutRedirect(bool enable) {
        if (enable) {
            old = std::cout.rdbuf();
            std::cout.rdbuf(nullptr);
        }
    }
    ~CoutRedirect() {
        if (old)
            std::cout.rdbuf(old);
    }
};

// Convert a tetwild mesh to numpy vertex and int arrays
nb::tuple extractMeshDataNumpy(const floatTetWild::Mesh &mesh, bool vtk_ordering) {
    std::map<int, int> old2new;
    int nv = 0;
    for (int i = 0; i < mesh.tet_vertices.size(); ++i) {
        if (!mesh.tet_vertices[i].is_removed) {
            old2new[i] = nv++;
        }
    }

    int nt = 0;
    for (const auto &t : mesh.tets)
        if (!t.is_removed)
            ++nt;

    // make numpy arrays and accessors
    NDArray<double, 2> V = MakeNDArray<double, 2>({nv, 3});
    NDArray<int, 2> T = MakeNDArray<int, 2>({nt, 4});
    double *v = V.data();
    int *t = T.data();

    int vi = 0;
    for (const auto &vert : mesh.tet_vertices) {
        if (vert.is_removed)
            continue;
        v[vi * 3 + 0] = vert.pos.x();
        v[vi * 3 + 1] = vert.pos.y();
        v[vi * 3 + 2] = vert.pos.z();
        ++vi;
    }

    // VTK ordering is different than tetwild's ordering
    int ti = 0;
    if (vtk_ordering) {
        for (const auto &tet : mesh.tets) {
            if (tet.is_removed)
                continue;
            t[ti * 4 + 0] = old2new[tet.indices[0]];
            t[ti * 4 + 1] = old2new[tet.indices[1]];
            t[ti * 4 + 2] = old2new[tet.indices[3]];
            t[ti * 4 + 3] = old2new[tet.indices[2]];
            ++ti;
        }
    } else {
        for (const auto &tet : mesh.tets) {
            if (tet.is_removed)
                continue;
            t[ti * 4 + 0] = old2new[tet.indices[0]];
            t[ti * 4 + 1] = old2new[tet.indices[1]];
            t[ti * 4 + 2] = old2new[tet.indices[2]];
            t[ti * 4 + 3] = old2new[tet.indices[3]];
            ++ti;
        }
    }

    return nb::make_tuple(V, T);
}

// convert a numpy array to a geogram vector
template <typename T> GEO::vector<T> array_to_geo_vector(NDArray<T, 2> array) {
    int sz = array.shape(0) * 3;
    GEO::vector<T> geo_vec(sz);

    // numpy arrays are contiguous
    std::memcpy(geo_vec.data(), array.data(), sz * sizeof(T));

    return geo_vec;
}

nb::tuple Tetrahedralize(
    NDArray<double, 2> vertices_arr,
    NDArray<unsigned int, 2> faces_arr,
    bool optimize,
    bool skip_simplify,
    double edge_length_r,
    double edge_length_abs,
    double epsilon,
    double stop_energy,
    bool coarsen,
    int num_threads,
    int num_opt_iter,
    int loglevel,
    bool quiet,
    bool vtk_ordering,
    bool disable_filtering,
    NDArray<double, 2> bg_vertices,
    NDArray<int, 2> bg_tets,
    NDArray<double, 1> bg_values) {
    using namespace floatTetWild;
    using namespace Eigen;
    GEO::initialize();
    CoutRedirect cout_guard(quiet);

    if (!quiet) {
        std::cout << "Starting tetrahedralization..." << std::endl;
    }

    // Initialize placeholders for flags and epsr_flags, if needed
    std::vector<int> flags;

    // Initialize GEO::Mesh and load the mesh data into it
    GEO::Mesh sf_mesh;
    if (!quiet) {
        std::cout << "Loading mesh..." << std::endl;
    }

    // geo_index_t is unsigned int
    GEO::vector<double> vertices_vec = array_to_geo_vector(vertices_arr);
    GEO::vector<geo_index_t> faces_vec = array_to_geo_vector(faces_arr);
    sf_mesh.facets.assign_triangle_mesh(3, vertices_vec, faces_vec, false);
    if (sf_mesh.cells.nb() != 0 && sf_mesh.facets.nb() == 0) {
        sf_mesh.cells.compute_borders();
    }
    if (!quiet) {
        sf_mesh.show_stats("I/O");
        std::cout << "Loaded mesh data into GEO::Mesh." << std::endl;
    }

    // Initialize AABBWrapper with the loaded GEO::Mesh for collision checking
    AABBWrapper tree(sf_mesh);
    if (!quiet) {
        std::cout << "Initialized AABBWrapper." << std::endl;
    }

    // Create an instance of Mesh to hold the output tetrahedral mesh
    Mesh mesh;
    if (!quiet) {
        std::cout << "Created Mesh instance for output." << std::endl;
    }

    // Prepare a vector to track the insertion status of faces
    std::vector<Eigen::Matrix<double, 3, 1>> input_points;
    std::vector<Eigen::Matrix<int, 3, 1>> input_faces;

    input_points.resize(sf_mesh.vertices.nb());
    for (size_t i = 0; i < input_points.size(); i++)
        input_points[i] << (sf_mesh.vertices.point(i))[0], (sf_mesh.vertices.point(i))[1],
            (sf_mesh.vertices.point(i))[2];

    input_faces.resize(sf_mesh.facets.nb());
    for (size_t i = 0; i < input_faces.size(); i++)
        input_faces[i] << sf_mesh.facets.vertex(i, 0), sf_mesh.facets.vertex(i, 1),
            sf_mesh.facets.vertex(i, 2);

    Parameters &params = mesh.params;
    params.eps_rel = epsilon;
    params.ideal_edge_length_rel = edge_length_r;
    params.ideal_edge_length_abs = edge_length_abs;
    params.stop_energy = stop_energy;
    params.coarsen = coarsen;
    params.is_quiet = quiet;
    params.max_its = num_opt_iter;
    params.disable_filtering = disable_filtering;

    if (!params.init(tree.get_sf_diag())) {
        throw std::runtime_error("FTetWildWrapper.cpp: Parameters initialization failed");
    }

    // Sizing field, from a background tet mesh carrying a target edge length
    // per vertex. fTetWild declares params.get_sizing_field_value and calls it
    // from apply_sizingfield(), but never assigns it anywhere -- the lookup
    // that main.cpp's --bg-mesh implies is commented out in
    // MeshImprovement.cpp, so setting apply_sizing_field alone would call an
    // empty std::function. Supply it here, following that original: find the
    // background tet containing the point and interpolate the target length
    // over its four vertices. These locals outlive the optimization() call
    // below, which is the only thing that reads the field.
    GEO::Mesh bg_mesh;
    std::unique_ptr<GEO::MeshCellsAABB> bg_tree;
    Eigen::VectorXd bg_V;
    Eigen::VectorXi bg_T;
    Eigen::VectorXd bg_vals;
    const size_t n_bg_verts = bg_vertices.shape(0);
    const size_t n_bg_tets = bg_tets.shape(0);
    if (n_bg_verts > 0 && n_bg_tets > 0) {
        bg_V.resize(static_cast<Eigen::Index>(n_bg_verts * 3));
        std::memcpy(bg_V.data(), bg_vertices.data(), n_bg_verts * 3 * sizeof(double));
        bg_vals.resize(static_cast<Eigen::Index>(n_bg_verts));
        std::memcpy(bg_vals.data(), bg_values.data(), n_bg_verts * sizeof(double));

        bg_T.resize(static_cast<Eigen::Index>(n_bg_tets * 4));
        const int *tet_src = bg_tets.data();
        for (size_t i = 0; i < n_bg_tets * 4; ++i) {
            if (tet_src[i] < 0 || static_cast<size_t>(tet_src[i]) >= n_bg_verts) {
                throw std::runtime_error(
                    "FTetWildWrapper.cpp: sizing field tetrahedron index out of range");
            }
            bg_T(static_cast<Eigen::Index>(i)) = tet_src[i];
        }

        bg_mesh.vertices.create_vertices(static_cast<int>(n_bg_verts));
        for (size_t i = 0; i < n_bg_verts; ++i) {
            GEO::vec3 &p = bg_mesh.vertices.point(static_cast<GEO::index_t>(i));
            for (int j = 0; j < 3; ++j)
                p[j] = bg_V(static_cast<Eigen::Index>(i * 3 + j));
        }
        bg_mesh.cells.create_tets(static_cast<int>(n_bg_tets));
        for (size_t i = 0; i < n_bg_tets; ++i) {
            for (int j = 0; j < 4; ++j) {
                bg_mesh.cells.set_vertex(
                    static_cast<GEO::index_t>(i),
                    static_cast<GEO::index_t>(j),
                    static_cast<GEO::index_t>(bg_T(static_cast<Eigen::Index>(i * 4 + j))));
            }
        }

        // The const overload builds the tree indirectly, leaving the mesh
        // ordering alone, so containing_tet() still returns an index into
        // bg_T rather than into a permuted copy.
        const GEO::Mesh &bg_mesh_ref = bg_mesh;
        bg_tree.reset(new GEO::MeshCellsAABB(bg_mesh_ref));

        // Only these two matter. The V/T/values_sizing_field members that
        // main.cpp fills in are read by nothing in fTetWild, so setting them
        // here would be dead weight.
        params.apply_sizing_field = true;
        params.get_sizing_field_value =
            [&bg_tree, &bg_V, &bg_T, &bg_vals](const floatTetWild::Vector3 &p) -> double {
            const GEO::index_t t_id = bg_tree->containing_tet(GEO::vec3(p[0], p[1], p[2]));
            if (t_id == GEO::MeshCellsAABB::NO_TET) {
                // Outside the background mesh. Zero leaves the point at the
                // global ideal edge length rather than collapsing it.
                return 0.0;
            }

            const Eigen::Index base = static_cast<Eigen::Index>(t_id) * 4;
            std::array<floatTetWild::Vector3, 4> vs;
            for (int j = 0; j < 4; ++j) {
                const Eigen::Index v = bg_T(base + j) * 3;
                vs[j] = floatTetWild::Vector3(bg_V(v), bg_V(v + 1), bg_V(v + 2));
            }

            // Barycentric weight of each vertex is the point's distance from
            // the opposite face over that vertex's distance from it.
            double value = 0.0;
            for (int j = 0; j < 4; ++j) {
                const floatTetWild::Vector3 n =
                    ((vs[(j + 1) % 4] - vs[j]).cross(vs[(j + 2) % 4] - vs[j])).normalized();
                const double d = (vs[(j + 3) % 4] - vs[j]).dot(n);
                if (d == 0)
                    continue;
                value +=
                    std::fabs((p - vs[j]).dot(n) / d) * bg_vals(bg_T(base + (j + 3) % 4));
            }
            return value;
        };
    }

    floatTetWild::Logger::init(!params.is_quiet, params.log_path);
    params.log_level = loglevel;
    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    // Set up threading
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }
    // IGL has issues with a nested for loop and oversubscription, see
    // https://github.com/libigl/libigl/issues/2412
    igl::default_num_threads(std::ceil(std::sqrt(num_threads)));
    params.num_threads = num_threads;
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    if (!quiet) {
        std::cout << "TBB threads " << num_threads << std::endl;
    }
    tbb::global_control parallelism_limit(
        tbb::global_control::max_allowed_parallelism, num_threads);
    tbb::global_control stack_size_limit(tbb::global_control::thread_stack_size, stack_size);

    std::vector<int> input_tags;
    if (input_tags.size() != input_faces.size()) {
        input_tags.resize(input_faces.size());
        std::fill(input_tags.begin(), input_tags.end(), 0);
    }
    // bool skip_simplify = false;
    simplify(input_points, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_points, input_faces, mesh);

    // Perform tetrahedralization
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    if (!quiet) {
        std::cout << "Starting tetrahedralization..." << std::endl;
    }
    FloatTetDelaunay::tetrahedralize(input_points, input_faces, tree, mesh, is_face_inserted);
    if (!quiet) {
        std::cout << "Tetrahedralization performed." << std::endl;
    }

    insert_triangles(
        input_points, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    if (optimize) {
        optimization(
            input_points,
            input_faces,
            input_tags,
            is_face_inserted,
            mesh,
            tree,
            {{1, 1, 1, 1}});
    }
    correct_tracked_surface_orientation(mesh, tree);

    // filter elements
    if (params.smooth_open_boundary) {
        smooth_open_boundary(mesh, tree);
        for (auto &t : mesh.tets) {
            if (t.is_outside)
                t.is_removed = true;
        }
    } else {
        if (!params.disable_filtering) {
            if (params.use_floodfill) {
                filter_outside_floodfill(mesh);
            } else if (params.use_input_for_wn) {
                filter_outside(mesh, input_points, input_faces);
            } else
                filter_outside(mesh);
        }
    }

    if (!quiet) {
        std::cout << "Tetrahedralization completed. Extracting mesh data..." << std::endl;
    }

    return extractMeshDataNumpy(mesh, vtk_ordering);
}

nb::tuple TetrahedralizeCSG(
    const std::string &csg_file,
    float epsilon,
    float edge_length_r,
    float stop_energy,
    bool coarsen,
    int num_threads,
    int log_level,
    bool vtk_ordering) {
    GEO::initialize();

    // Initialize mesh and parameters
    floatTetWild::Mesh mesh;
    floatTetWild::Parameters &params = mesh.params;
    params.eps_rel = epsilon;
    params.ideal_edge_length_rel = edge_length_r;
    params.stop_energy = stop_energy;
    params.coarsen = coarsen;
    // params.use_general_wn = true;

    // Set up threading
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }
    // IGL has issues with a nested for loop and oversubscription, see
    // https://github.com/libigl/libigl/issues/2412
    igl::default_num_threads(std::ceil(std::sqrt(num_threads)));
    params.num_threads = num_threads;
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    std::cout << "TBB threads " << num_threads << std::endl;
    tbb::global_control parallelism_limit(
        tbb::global_control::max_allowed_parallelism, num_threads);
    tbb::global_control stack_size_limit(tbb::global_control::thread_stack_size, stack_size);

    floatTetWild::Logger::init(!params.is_quiet, params.log_path);
    params.log_level = log_level;
    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    // Load CSG tree
    json csg_tree = json({});
    std::ifstream file(csg_file);
    if (file.is_open()) {
        file >> csg_tree;
        file.close();
    } else {
        throw std::runtime_error("Unable to open CSG file: " + csg_file);
    }

    // Load and merge meshes from CSG tree
    json tree_with_ids;
    std::vector<std::string> meshes;
    CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);

    // Pre-check all input meshes for non-finite (NaN/Inf) coordinates.
    nb::print("Pre-checking CSG input meshes for NaN/Inf values...");
    for (const auto &mesh_file : meshes) {
        std::vector<Vector3> temp_vertices;
        std::vector<Vector3i> temp_faces;
        std::vector<int> temp_tags;
        GEO::Mesh temp_sf_mesh; // Dummy, but load_and_merge requires it

        // We use load_and_merge with a single-item list
        if (!CSGTreeParser::load_and_merge(
                {mesh_file}, // Load just this one file
                temp_vertices,
                temp_faces,
                temp_sf_mesh,
                temp_tags)) {
            throw std::runtime_error("Failed to pre-load mesh for checking: " + mesh_file);
        }

        // The check for NaN/Inf
        for (const auto &v : temp_vertices) {
            if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
                throw std::runtime_error(
                    "FATAL: Input mesh file contains non-finite (NaN/Inf) vertex "
                    "coordinates: " +
                    mesh_file);
            }
        }
    }
    nb::print("All input meshes passed check.");

    std::vector<Vector3> input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int> input_tags;
    GEO::Mesh sf_mesh;
    if (!CSGTreeParser::load_and_merge(
            meshes, input_vertices, input_faces, sf_mesh, input_tags))
        throw std::runtime_error("Failed to load and merge meshes from CSG tree");

    // Initialize AABBWrapper
    AABBWrapper tree(sf_mesh);
    if (!params.init(tree.get_sf_diag())) {
        throw std::runtime_error("Parameters initialization failed");
    }

    // Preprocessing
    bool skip_simplify = false;
    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);

    // Tetrahedralization
    std::vector<bool> is_face_inserted(input_faces.size(), false);
    FloatTetDelaunay::tetrahedralize(
        input_vertices, input_faces, tree, mesh, is_face_inserted);

    // Insert triangles and optimize
    insert_triangles(
        input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
    optimization(
        input_vertices,
        input_faces,
        input_tags,
        is_face_inserted,
        mesh,
        tree,
        {{1, 1, 1, 1}});

    std::cout << "optimization finished" << std::endl;
    // Correct surface orientation
    correct_tracked_surface_orientation(mesh, tree);
    std::cout << "correct surface orientation finished" << std::endl;

    // Apply Boolean operations
    boolean_operation(mesh, tree_with_ids, meshes);
    std::cout << "boolean operation finished " << std::endl;

    // Extract data
    nb::tuple result = extractMeshDataNumpy(mesh, vtk_ordering);
    NDArray<double, 2> np_vertices = nb::cast<NDArray<double, 2>>(result[0]);
    NDArray<int, 2> np_cells = nb::cast<NDArray<int, 2>>(result[1]);

    // Extract marker
    std::vector<int> marker;
    for (const auto &tet : mesh.tets) {
        if (!tet.is_removed) {
            marker.push_back(tet.scalar);
        }
    }

    // Extract markers
    int num_cells = np_cells.shape(0);
    NDArray<int, 1> np_markers = MakeNDArray<int, 1>({num_cells});
    int *m_data = np_markers.data();
    for (size_t i = 0; i < num_cells; ++i) {
        m_data[i] = marker[i];
    }

    nb::print("Tetrahedralization process from CSG completed successfully.");
    return nb::make_tuple(np_vertices, np_cells, np_markers);
}

NB_MODULE(PyfTetWildWrapper, m) {
    m.doc() = "Nanobind wrapper for FloatTetWild mesh tetrahedralization";
    m.def(
        "tetrahedralize_mesh",
        &Tetrahedralize,
        R"doc(
Tetrahedralizes a mesh given vertices and faces arrays, returning numpy arrays
of tetrahedra and points.)doc");

    m.def(
        "tetrahedralize_csg",
        &TetrahedralizeCSG,
        R"doc(
Tetrahedralizes a CSG tree from a JSON file, returning numpy arrays of
vertices, cells, and markers.)doc",
        nb::arg("csg_file"),
        nb::arg("epsilon"),
        nb::arg("edge_length_r"),
        nb::arg("stop_energy"),
        nb::arg("coarsen"),
        nb::arg("num_threads"),
        nb::arg("log_level"),
        nb::arg("vtk_ordering"));
}
