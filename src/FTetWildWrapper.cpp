#include <floattetwild/Mesh.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/MeshIO.hpp>

#include <Eigen/Dense>
#include <vector>

// You must include other necessary headers and namespaces used in your project

void tetrahedralizeAndWriteMesh(const std::vector<Eigen::Vector3d>& vertices,
                                const std::vector<Eigen::Vector3i>& faces,
                                const std::string& outputPath) {
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

    // Write the tetrahedralized mesh to a file
    std::vector<Scalar> colors; // Optional: For visualizing quality or other metrics
    MeshIO::write_mesh(outputPath, mesh, false, colors);
}

// Example usage:
int main() {
    std::vector<Eigen::Vector3d> vertices = {/* Fill with your vertices */};
    std::vector<Eigen::Vector3i> faces = {/* Fill with your faces */};
    std::string outputPath = "output_mesh.msh";

    tetrahedralizeAndWriteMesh(vertices, faces, outputPath);

    return 0;
}
