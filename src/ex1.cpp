#include "include/ex1.h"

// ======================================================================
// PROBLEM 1
// ========================================================================
void Exercise1::solve_sparse_linear_system()
{
    SMatd A(5, 5);
    std::vector<T> tripletList;
    VecXd b(5);

    // TODO: ADD YOUR CODE HERE

    // Adapted from https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
    tripletList.push_back(T(0, 1, 1));
    tripletList.push_back(T(0, 3, -2));
    tripletList.push_back(T(1, 0, -1));
    tripletList.push_back(T(1, 2, 3));
    tripletList.push_back(T(1, 4, 4));
    tripletList.push_back(T(2, 3, 5));
    tripletList.push_back(T(2, 4, 2));
    tripletList.push_back(T(3, 0, -1));
    tripletList.push_back(T(3, 1, 3));
    tripletList.push_back(T(4, 2, 1));
    tripletList.push_back(T(4, 3, 1));
    tripletList.push_back(T(4, 4, 1));

    A.setFromTriplets(tripletList.begin(), tripletList.end());

    b << 1., 10., 8., 0., 3.;

    A.makeCompressed();
    VecXd x(5);
    Eigen::SparseLU<SMatd, Eigen::COLAMDOrdering<int>> solver;
    // fill A and b;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(A);
    // Compute the numerical factorization
    solver.factorize(A);
    // Use the factors to solve the linear system
    x = solver.solve(b);

    std::cout << x << std::endl;
}

// ======================================================================
// PROBLEM 2
// ========================================================================
void Exercise1::generate_plot(int function_id)
{

    switch (function_id)
    {
    case 0:
    {

        // TODO: ADD YOUR CODE FOR 2(a) HERE
        // geometry.requireEdgeLengths()
        geometry->requireVertexPositions();
        surface::VertexData<double> vHeights(*gc_mesh);

        for (Vertex v : gc_mesh->vertices())
        {
            vHeights[v] = geometry->vertexPositions[v][1];
        }

        auto sc = ps_mesh->addVertexScalarQuantity("Vertex Heights", vHeights);
        sc->setEnabled(true);
        break;
    }
    case 1:
    {

        // TODO: ADD YOUR CODE FOR 2(b) HERE
        geometry->requireEdgeLengths();
        surface::EdgeData<double> edgeLengths(*gc_mesh);

        for (Edge e : gc_mesh->edges())
        {
            edgeLengths[e] = geometry->edgeLengths[e];
        }

        auto sc = ps_mesh->addEdgeScalarQuantity("Edge Lengths", edgeLengths);
        sc->setEnabled(true);
        break;
    }
    case 2:
    {

        // TODO: ADD YOUR CODE FOR 2(c) HERE
        geometry->requireFaceNormals();
        surface::FaceData<Vector3> FaceNormals(*gc_mesh);

        for (Face f : gc_mesh->faces())
        {
            FaceNormals[f] = geometry->faceNormals[f];
            FaceNormals[f][0] = (FaceNormals[f][0] + 1) / 2;
            FaceNormals[f][1] = (FaceNormals[f][1] + 1) / 2;
            FaceNormals[f][2] = (FaceNormals[f][2] + 1) / 2;
        }

        auto sc = ps_mesh->addFaceColorQuantity("Face Normals", FaceNormals);
        sc->setEnabled(true);
        break;
    }
    }
}
