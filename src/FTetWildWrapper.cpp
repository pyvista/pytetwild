#include <floattetwild/Mesh.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/MeshIO.hpp>

#include <Eigen/Dense>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <vector>
#include <floattetwild/Mesh.hpp>

std::pair<std::vector<std::array<float, 3>>, std::vector<std::array<int, 4>>> extractMeshData(const floatTetWild::Mesh& mesh) {
    std::vector<std::array<float, 3>> vertices;
    std::vector<std::array<int, 4>> tetrahedra;

    // Extract vertices
    for (const auto& vertex : mesh.tet_vertices) {
        if (!vertex.is_removed) {
            vertices.push_back({{static_cast<float>(vertex.pos.x()), static_cast<float>(vertex.pos.y()), static_cast<float>(vertex.pos.z())}});
        }
    }

    // Extract tetrahedra
    for (const auto& tet : mesh.tets) {
        if (!tet.is_removed) {
            tetrahedra.push_back({{tet.indices[0], tet.indices[1], tet.indices[2], tet.indices[3]}});
        }
    }

    return {vertices, tetrahedra};
}


std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector4i>> tetrahedralizeAndWriteMesh(
  const std::vector<Eigen::Vector3d>& vertices,
  const std::vector<Eigen::Vector3i>& faces) {
    using namespace floatTetWild;
    using namespace Eigen;

    // Convert input vertices and faces to the format expected by FloatTetWild
    std::vector<Vector3> input_vertices(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        input_vertices[i] = vertices[i].cast<Scalar>();
    }

    std::vector<Vector3i> input_faces(faces.begin(), faces.end());

    // Initialize placeholders for flags and epsr_flags, if needed
    std::vector<int> flags;
    std::vector<double> epsr_flags;

    // Initialize GEO::Mesh and load the mesh data into it
    GEO::Mesh geo_mesh;
    MeshIO::load_mesh(input_vertices, input_faces, geo_mesh, flags, epsr_flags);

    // Initialize AABBWrapper with the loaded GEO::Mesh for collision checking
    AABBWrapper tree(geo_mesh);

    // Create an instance of Mesh to hold the output tetrahedral mesh
    Mesh mesh;

    // Prepare a vector to track the insertion status of faces
    std::vector<bool> is_face_inserted(input_faces.size(), false);

    // Perform tetrahedralization
    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);

    return extractMeshData(mesh);
}

PYBIND11_MODULE(_wrapper, m) {
    m.doc() = "Pybind11 plugin for FloatTetWild mesh tetrahedralization";

    m.def("tetrahedralize_mesh", [](py::array_t<double> vertices, py::array_t<int> faces) {
        // Convert numpy arrays to Eigen matrices for vertices and faces
        auto verts = vertices.unchecked<2>(); // Access numpy array data without bounds checking
        auto fcs = faces.unchecked<2>();

        std::vector<Eigen::Vector3d> vert_vector(verts.shape(0));
        std::vector<Eigen::Vector3i> face_vector(fcs.shape(0));

        for (py::ssize_t i = 0; i < verts.shape(0); ++i) {
            vert_vector[i] = Eigen::Vector3d(verts(i, 0), verts(i, 1), verts(i, 2));
        }

        for (py::ssize_t i = 0; i < fcs.shape(0); ++i) {
            face_vector[i] = Eigen::Vector3i(fcs(i, 0), fcs(i, 1), fcs(i, 2));
        }

        // Call the tetrahedralization function
        auto [vertices_result, tetrahedra_result] = tetrahedralizeAndWriteMesh(vert_vector, face_vector);

        // Convert results back to numpy arrays
        size_t num_vertices = vertices_result.size();
        size_t num_tetrahedra = tetrahedra_result.size();

        // Prepare numpy array (points)
        py::array_t<float> np_vertices({num_vertices, 3});
        auto np_verts_access = np_vertices.mutable_unchecked<2>();
        for (size_t i = 0; i < num_vertices; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                np_verts_access(i, j) = vertices_result[i][j];
            }
        }

        // Prepare numpy array (tetrahedra)
        py::array_t<int> np_tetrahedra({num_tetrahedra, 4});
        auto np_tets_access = np_tetrahedra.mutable_unchecked<2>();
        for (size_t i = 0; i < num_tetrahedra; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                np_tets_access(i, j) = tetrahedra_result[i][j];
            }
        }

        return std::make_pair(np_vertices, np_tetrahedra);
    }, "Tetrahedralizes a mesh given vertices and faces arrays, returning numpy arrays of tetrahedra and points.");
}


