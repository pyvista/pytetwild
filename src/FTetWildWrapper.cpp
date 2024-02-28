#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/FloatTetwild.h>
#include <Eigen/Dense>

// Function to generate tetrahedral mesh using fTetWild
// vertices: Input array of surface vertices, each vertex is represented by 3 consecutive doubles (x, y, z)
// numVertices: Number of vertices in the input surface mesh
// faces: Input array of surface triangle faces, each face is represented by 3 consecutive integers (indices to vertices)
// numFaces: Number of faces in the input surface mesh
// tetPoints: Output array of tetrahedral mesh vertices, each vertex is represented by 3 consecutive doubles (x, y, z)
// numTetPoints: Output number of vertices in the tetrahedral mesh
// tetCells: Output array of tetrahedral mesh cells, each cell is represented by 4 consecutive integers (indices to tetPoints)
// numTetCells: Output number of cells in the tetrahedral mesh
void generateTetMesh(const double* vertices, int numVertices, const int* faces, int numFaces, double** tetPoints, int* numTetPoints, int** tetCells, int* numTetCells) {
    // Convert input arrays to Eigen matrices for fTetWild
    Eigen::MatrixXd V(numVertices, 3);
    Eigen::MatrixXi F(numFaces, 3);
    for(int i = 0; i < numVertices; ++i) {
        V(i, 0) = vertices[i * 3];
        V(i, 1) = vertices[i * 3 + 1];
        V(i, 2) = vertices[i * 3 + 2];
    }
    for(int i = 0; i < numFaces; ++i) {
        F(i, 0) = faces[i * 3];
        F(i, 1) = faces[i * 3 + 1];
        F(i, 2) = faces[i * 3 + 2];
    }

    // Parameters for tetrahedralization
    floatTetWild::Parameters params;
    // Set any required parameters here
    // For example, params.ideal_edge_length = 1.0;

    // Mesh object to store the output tetrahedral mesh
    floatTetWild::Mesh mesh;

    // Generate the tetrahedral mesh
    floatTetWild::tetrahedralize(V, F, mesh, params);

    // Extract tetrahedral mesh vertices and cells
    auto& vertices_out = mesh.get_vertices();
    auto& tets_out = mesh.get_tets();

    // Assuming memory for tetPoints and tetCells is allocated externally
    *numTetPoints = vertices_out.size();
    *numTetCells = tets_out.size();

    for(int i = 0; i < vertices_out.size(); ++i) {
        (*tetPoints)[i * 3] = vertices_out[i][0];
        (*tetPoints)[i * 3 + 1] = vertices_out[i][1];
        (*tetPoints)[i * 3 + 2] = vertices_out[i][2];
    }

    for(int i = 0; i < tets_out.size(); ++i) {
        if(tets_out[i].is_removed) continue; // Skip removed tets
        (*tetCells)[i * 4] = tets_out[i][0];
        (*tetCells)[i * 4 + 1] = tets_out[i][1];
        (*tetCells)[i * 4 + 2] = tets_out[i][2];
        (*tetCells)[i * 4 + 3] = tets_out[i][3];
    }
}
