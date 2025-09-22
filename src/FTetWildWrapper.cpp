#include <floattetwild/AABBWrapper.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/Types.hpp>
#include <thread>

#include <Eigen/Dense>
#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
namespace py = pybind11;

#include <cstdint>
#include <vector>

#include <geogram/api/defs.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>

using floatTetWild::AABBWrapper;
using floatTetWild::boolean_operation;
using floatTetWild::CSGTreeParser;
using floatTetWild::FloatTetDelaunay;
using floatTetWild::json;
using floatTetWild::Vector3;
using floatTetWild::Vector3i;

std::pair<std::vector<std::array<float, 3>>, std::vector<std::array<int, 4>>>
extractMeshData(const floatTetWild::Mesh &mesh) {
  std::vector<std::array<float, 3>> vertices;
  std::vector<std::array<int, 4>> tetrahedra;

  // Remap tet indices for removed vertices
  std::map<int, int> old_2_new;
  int cnt_v = 0;
  const auto skip_vertex = [&mesh](const int i) {
    return mesh.tet_vertices[i].is_removed;
  };
  for (int i = 0; i < mesh.tet_vertices.size(); i++) {
    if (!skip_vertex(i)) {
      old_2_new[i] = cnt_v;
      cnt_v++;
    }
  }

  // Extract vertices
  for (const auto &vertex : mesh.tet_vertices) {
    if (!vertex.is_removed) {
      vertices.push_back({{static_cast<float>(vertex.pos.x()),
                           static_cast<float>(vertex.pos.y()),
                           static_cast<float>(vertex.pos.z())}});
    }
  }

  // Extract tetrahedra
  for (const auto &tet : mesh.tets) {
    if (!tet.is_removed) {
      tetrahedra.push_back(
          {{old_2_new[tet.indices[0]], old_2_new[tet.indices[1]],
            old_2_new[tet.indices[2]], old_2_new[tet.indices[3]]}});
    }
  }

  return {vertices, tetrahedra};
}

std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<int, 4>>,
           std::vector<int>>
extractMeshDataMarker(const floatTetWild::Mesh &mesh) {
  std::vector<std::array<float, 3>> vertices;
  std::vector<std::array<int, 4>> tetrahedra;
  std::vector<int> marker;

  // Remap tet indices for removed vertices
  std::map<int, int> old_2_new;
  int cnt_v = 0;
  const auto skip_vertex = [&mesh](const int i) {
    return mesh.tet_vertices[i].is_removed;
  };
  for (int i = 0; i < mesh.tet_vertices.size(); i++) {
    if (!skip_vertex(i)) {
      old_2_new[i] = cnt_v;
      cnt_v++;
    }
  }

  // Extract vertices
  for (const auto &vertex : mesh.tet_vertices) {
    if (!vertex.is_removed) {
      vertices.push_back({{static_cast<float>(vertex.pos.x()),
                           static_cast<float>(vertex.pos.y()),
                           static_cast<float>(vertex.pos.z())}});
    }
  }

  // Extract tetrahedra
  for (const auto &tet : mesh.tets) {
    if (!tet.is_removed) {
      tetrahedra.push_back(
          {{old_2_new[tet.indices[0]], old_2_new[tet.indices[1]],
            old_2_new[tet.indices[2]], old_2_new[tet.indices[3]]}});
    }
  }
  // Extract marker
  for (const auto &tet : mesh.tets) {
    if (!tet.is_removed) {
      marker.push_back(old_2_new[tet.scalar]);
    }
  }

  return {vertices, tetrahedra, marker};
}

template <typename T> GEO::vector<T> array_to_geo_vector(py::array_t<T> array) {
  auto info = array.request();
  T *ptr = static_cast<T *>(info.ptr);
  GEO::vector<T> geo_vec(info.size); // Allocate space for info.size elements.

  // Populate the GEO::vector<T> with elements from the array.
  for (size_t i = 0; i < info.size; ++i) {
    geo_vec.data()[i] = ptr[i]; // Copy each element.
  }

  return geo_vec;
}

std::pair<std::vector<std::array<float, 3>>, std::vector<std::array<int, 4>>>
tetrahedralize(GEO::vector<double> &vertices, GEO::vector<geo_index_t> &faces,
               bool optimize, float scale_fac) {
  using namespace floatTetWild;
  using namespace Eigen;

  std::cout << "Starting tetrahedralization..." << std::endl;

  // Initialize placeholders for flags and epsr_flags, if needed
  std::vector<int> flags;

  // Initialize GEO::Mesh and load the mesh data into it
  GEO::Mesh sf_mesh;
  std::cout << "Loading mesh..." << std::endl;

  sf_mesh.facets.assign_triangle_mesh(3, vertices, faces, false);
  if (sf_mesh.cells.nb() != 0 && sf_mesh.facets.nb() == 0) {
    sf_mesh.cells.compute_borders();
  }
  sf_mesh.show_stats("I/O");

  // std::cout << "Vertices (Points):" << std::endl;
  // for (size_t i = 0; i < sf_mesh.vertices.nb(); ++i) {
  //   const GEO::vec3 &p =
  //       sf_mesh.vertices.point(i); // Access the point at index i
  //   std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " <<
  //   p[2]
  //             << ")" << std::endl;
  // }

  // Print all facets (faces)
  // std::cout << "\nFacets (Faces):" << std::endl;
  // for (size_t i = 0; i < sf_mesh.facets.nb(); ++i) {
  //   std::cout << "Face " << i << ":";
  //   for (size_t j = 0; j < sf_mesh.facets.nb_vertices(i); ++j) {
  //     // Print each vertex index that makes up the facet
  //     std::cout << " " << sf_mesh.facets.vertex(i, j);
  //   }
  //   std::cout << std::endl;
  // }

  GEO::mesh_reorder(sf_mesh, GEO::MESH_ORDER_MORTON);
  std::cout << "Loaded mesh data into GEO::Mesh." << std::endl;

  // Initialize AABBWrapper with the loaded GEO::Mesh for collision checking
  AABBWrapper tree(sf_mesh);
  std::cout << "Initialized AABBWrapper." << std::endl;

  // Create an instance of Mesh to hold the output tetrahedral mesh
  Mesh mesh;
  std::cout << "Created Mesh instance for output." << std::endl;

  // Prepare a vector to track the insertion status of faces
  std::vector<Eigen::Matrix<double, 3, 1>> input_points;
  std::vector<Eigen::Matrix<int, 3, 1>> input_faces;

  input_points.resize(sf_mesh.vertices.nb());
  for (size_t i = 0; i < input_points.size(); i++)
    input_points[i] << (sf_mesh.vertices.point(i))[0],
        (sf_mesh.vertices.point(i))[1], (sf_mesh.vertices.point(i))[2];

  input_faces.resize(sf_mesh.facets.nb());
  for (size_t i = 0; i < input_faces.size(); i++)
    input_faces[i] << sf_mesh.facets.vertex(i, 0), sf_mesh.facets.vertex(i, 1),
        sf_mesh.facets.vertex(i, 2);

  Parameters &params = mesh.params;
  if (!params.init(tree.get_sf_diag() * scale_fac)) {
    throw std::runtime_error(
        "FTetWildWrapper.cpp: Parameters initialization failed");
  }
  const size_t MB = 1024 * 1024;
  const size_t stack_size = 64 * MB;
  unsigned int max_threads = std::numeric_limits<unsigned int>::max();
  unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
  num_threads = std::min(max_threads, num_threads);
  params.num_threads = num_threads;

  std::vector<int> input_tags;
  if (input_tags.size() != input_faces.size()) {
    input_tags.resize(input_faces.size());
    std::fill(input_tags.begin(), input_tags.end(), 0);
  }
  bool skip_simplify = false;
  simplify(input_points, input_faces, input_tags, tree, params, skip_simplify);
  tree.init_b_mesh_and_tree(input_points, input_faces, mesh);

  // Perform tetrahedralization
  std::vector<bool> is_face_inserted(input_faces.size(), false);
  std::cout << "Starting tetrahedralization..." << std::endl;
  FloatTetDelaunay::tetrahedralize(input_points, input_faces, tree, mesh,
                                   is_face_inserted);
  std::cout << "Tetrahedralization performed." << std::endl;

  insert_triangles(input_points, input_faces, input_tags, mesh,
                   is_face_inserted, tree, false);
  if (optimize) {
    optimization(input_points, input_faces, input_tags, is_face_inserted, mesh,
                 tree, {{1, 1, 1, 1}});
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

  std::cout << "Tetrahedralization completed. Extracting mesh data..."
            << std::endl;

  return extractMeshData(mesh);
}

PYBIND11_MODULE(PyfTetWildWrapper, m) {
  m.doc() = "Pybind11 plugin for FloatTetWild mesh tetrahedralization";

  m.def(
      "tetrahedralize_mesh",
      [](py::array_t<double> vertices, py::array_t<unsigned int> faces,
         bool optimize, float scale_fac) {
        // GEO::Logger* geo_logger = GEO::Logger::instance();
        // geo_logger-->initialize();
        GEO::initialize();

        scale_fac *= 20; // already treated as 1/20 internally

        py::print("Starting tetrahedralization process...");

        // Convert numpy arrays to vectors and call the tetrahedralization
        // function
        GEO::vector<double> vertices_vec = array_to_geo_vector(vertices);
        GEO::vector<geo_index_t> faces_vec = array_to_geo_vector(faces);
        auto result =
            tetrahedralize(vertices_vec, faces_vec, optimize, scale_fac);
        auto vertices_result = result.first;
        auto tetrahedra_result = result.second;
        py::print("Tetrahedralization complete.");

        // Convert results back to numpy arrays
        size_t num_vertices = vertices_result.size();
        size_t num_tetrahedra = tetrahedra_result.size();
        py::print("Number of vertices:", num_vertices);
        py::print("Number of tetrahedra:", num_tetrahedra);

        // Prepare numpy array (points)
        size_t shape[2]{num_vertices, 3};
        auto np_vertices = py::array_t<float>(shape);
        auto np_verts_access = np_vertices.mutable_unchecked<2>();
        for (size_t i = 0; i < num_vertices; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            np_verts_access(i, j) = vertices_result[i][j];
          }
        }
        py::print("Prepared numpy array for points.");

        // Prepare numpy array (tetrahedra)
        size_t shape_tet[2]{num_tetrahedra, 4};
        auto np_tetrahedra = py::array_t<int>(shape_tet);
        auto np_tets_access = np_tetrahedra.mutable_unchecked<2>();
        for (size_t i = 0; i < num_tetrahedra; ++i) {
          for (size_t j = 0; j < 4; ++j) {
            np_tets_access(i, j) = tetrahedra_result[i][j];
          }
        }
        py::print("Prepared numpy array for tetrahedra.");

        py::print("Tetrahedralization process completed successfully.");
        return std::make_pair(np_vertices, np_tetrahedra);
      },
      "Tetrahedralizes a mesh given vertices and faces arrays, returning numpy "
      "arrays of tetrahedra and points.");

  m.def(
      "tetrahedralize_csg",
      [](const std::string &csg_file, float epsilon, float edge_length_r,
         float stop_energy, bool coarsen) {
        GEO::initialize();

        // Initialize mesh and parameters
        floatTetWild::Mesh mesh;
        floatTetWild::Parameters &params = mesh.params;
        params.eps_rel = epsilon;
        params.ideal_edge_length_rel = edge_length_r;
        params.stop_energy = stop_energy;
        params.coarsen = coarsen;

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

        std::vector<Vector3> input_vertices;
        std::vector<Vector3i> input_faces;
        std::vector<int> input_tags;
        GEO::Mesh sf_mesh;
        if (!CSGTreeParser::load_and_merge(meshes, input_vertices, input_faces,
                                           sf_mesh, input_tags))
          throw std::runtime_error(
              "Failed to load and merge meshes from CSG tree");

        // Initialize AABBWrapper
        AABBWrapper tree(sf_mesh);
        if (!params.init(tree.get_sf_diag())) {
          throw std::runtime_error("Parameters initialization failed");
        }

        // Set up threading
        unsigned int max_threads = std::numeric_limits<unsigned int>::max();
        unsigned int num_threads =
            std::max(1u, std::thread::hardware_concurrency());
        num_threads = std::min(max_threads, num_threads);
        params.num_threads = num_threads;

        // Preprocessing
        bool skip_simplify = false;
        simplify(input_vertices, input_faces, input_tags, tree, params,
                 skip_simplify);
        tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);

        // Tetrahedralization
        std::vector<bool> is_face_inserted(input_faces.size(), false);
        FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree,
                                         mesh, is_face_inserted);

        // Insert triangles and optimize
        insert_triangles(input_vertices, input_faces, input_tags, mesh,
                         is_face_inserted, tree, false);
        optimization(input_vertices, input_faces, input_tags, is_face_inserted,
                     mesh, tree, {{1, 1, 1, 1}});

        // Correct surface orientation
        correct_tracked_surface_orientation(mesh, tree);

        // Apply Boolean operations
        boolean_operation(mesh, tree_with_ids, meshes);

        // Extract surface data
        auto result = extractMeshDataMarker(mesh);
        auto vertices = std::get<0>(result);
        auto cells = std::get<1>(result);
        auto marker = std::get<2>(result);

        // Convert to numpy arrays
        size_t num_vertices = vertices.size();
        size_t num_cells = cells.size();

        size_t v_shape[2] = {num_vertices, 3};
        auto np_vertices = py::array_t<float>(v_shape);
        auto v_access = np_vertices.mutable_unchecked<2>();
        for (size_t i = 0; i < num_vertices; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            v_access(i, j) = vertices[i][j];
          }
        }

        size_t c_shape[2] = {num_cells, 4};
        auto np_cells = py::array_t<int>(c_shape);
        auto c_access = np_cells.mutable_unchecked<2>();
        for (size_t i = 0; i < num_cells; ++i) {
          for (size_t j = 0; j < 4; ++j) {
            c_access(i, j) = cells[i][j];
          }
        }

        size_t m_shape[1] = {num_cells};
        auto np_markers = py::array_t<int>(m_shape);
        auto m_access = np_markers.mutable_unchecked<1>();
        for (size_t i = 0; i < num_cells; ++i) {
          m_access(i) = marker[i];
        }
        py::print(
            "Tetrahedralization process from CSG completed successfully.");
        return std::make_tuple(np_vertices, np_cells, np_markers);
      },
      "Tetrahedralizes a CSG tree from a JSON file, returning numpy arrays of "
      "vertices, cells, and markers.",
      py::arg("csg_file"), py::arg("epsilon"), py::arg("edge_length_r"),
      py::arg("stop_energy"), py::arg("coarsen"));
}
