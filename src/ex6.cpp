#include "include/ex6.h"
#include "geometrycentral/numerical/linear_algebra_types.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/utilities/vector3.h"
#include <iostream>

// ======================================================================
// EXERCISE 1.1
// ========================================================================
void Exercise6::uniform_smooth()
{

    // TODO: ADD YOUR CODE HERE
    std::vector<Vector3> update;
    geometry->requireEdgeCotanWeights();
    geometry->requireVertexPositions();
    for (Vertex v : gc_mesh->vertices())
    {
        Point p = geometry->vertexPositions[v];
        if (!v.isBoundary())
        {
            int ncount = 0;
            Point sum = Vector3{0., 0., 0.};
            for (Vertex vn : v.adjacentVertices())
            {
                ncount += 1;
                sum += (geometry->vertexPositions[vn] - p);
            }
            p = p + 0.5 * (sum / ncount);
            update.push_back({p[0], p[1], p[2]});
        }
        else
        {
            update.push_back({p[0], p[1], p[2]});
        }
    }
    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        geometry->vertexPositions[j] = update[j];
    }
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}

// ======================================================================
// EXERCISE 1.2
// ========================================================================
void Exercise6::cotan_laplacian_smooth()
{

    // TODO: ADD YOUR CODE HERE
    std::vector<Vector3> update;
    geometry->requireEdgeCotanWeights();
    geometry->requireVertexPositions();
    for (Vertex v : gc_mesh->vertices())
    {
        Point p = geometry->vertexPositions[v];
        if (!v.isBoundary())
        {
            double ncount = 0;
            Point sum = Vector3{0., 0., 0.};
            for (Edge e : v.adjacentEdges())
            {
                double w = geometry->edgeCotanWeights[e];
                ncount += w;
                sum += w * (geometry->vertexPositions[e.otherVertex(v)] - p);
            }
            p = p + 0.5 * (sum / ncount);
        }
        update.push_back({p[0], p[1], p[2]});
    }
    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        geometry->vertexPositions[j] = update[j];
    }
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}

// ======================================================================
// EXERCISE 2
// ========================================================================
void Exercise6::implicit_smooth()
{

    // TODO: ADD YOUR CODE HERE
    Eigen::SparseMatrix<double> Dinv;
    std::vector<T> DtripletList;
    Eigen::SparseMatrix<double> M;
    std::vector<T> MtripletList;
    Eigen::MatrixXd current(gc_mesh->nVertices(), 3);
    geometry->requireEdgeCotanWeights();
    geometry->requireVertexPositions();
    Eigen::VectorXd d;
    for (Vertex v : gc_mesh->vertices())
    {
        int i = geometry->vertexIndices[v];
        Point p = geometry->vertexPositions[v];
        Eigen::Vector3d r(p[0], p[1], p[2]);

        current(i, 0) = p[0];
        current(i, 1) = p[1];
        current(i, 2) = p[2];
        // current.row(i) = r;

        DtripletList.push_back(T(i, i, (2 * geometry->vertexDualAreas[v]))); // take inverse by using reciprocal
        double totalw = 0;
        for (Edge e : v.adjacentEdges()) // Bug was that I was checking for i == j in this loop
        {
            double w = geometry->edgeCotanWeights[e];
            totalw += geometry->edgeCotanWeights[e];
            int j = geometry->vertexIndices[e.otherVertex(v)];
            MtripletList.push_back(T(i, j, w));
        }
        MtripletList.push_back(T(i, i, -1 * totalw));
    }
    Dinv.resize(DtripletList.size(), DtripletList.size());
    Dinv.setFromTriplets(DtripletList.begin(), DtripletList.end());
    M.resize(DtripletList.size(), DtripletList.size()); // should be same size as D
    M.setFromTriplets(MtripletList.begin(), MtripletList.end());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
    cg.compute(Dinv - (timestep * (M))); // timestep?
    Eigen::MatrixXd update = cg.solve((Dinv * current));
    // Eigen::MatrixXd update = (Eigen::MatrixXd(Dinv - (timestep * (M)))).colPivHouseholderQr().solve((Eigen::MatrixXd(Dinv * current))); // TEST
    // std::cout << "DEBUG1 " << current << std::endl; // sparse solvers give NaN for y and z, but general solver only gives NaN for z
    // std::cout << "DEBUG2 " << Dinv << std::endl;
    // std::cout << "DEBUG3 " << M << std::endl;
    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        geometry->vertexPositions[j] = {update(j, 0), update(j, 1), update(j, 2)};
    }
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}

// ======================================================================
// EXERCISE 3.1
// ========================================================================
void Exercise6::uniform_laplacian_enhance_feature()
{

    // TODO: ADD YOUR CODE HERE
    std::vector<Vector3> update;
    geometry->requireEdgeCotanWeights();
    geometry->requireVertexPositions();
    for (Vertex v : gc_mesh->vertices())
    {
        Point p = geometry->vertexPositions[v];
        if (!v.isBoundary())
        {
            int ncount = 0;
            Point sum = Vector3{0., 0., 0.};
            for (Vertex vn : v.adjacentVertices())
            {
                ncount += 1;
                sum += (geometry->vertexPositions[vn] - p);
            }
            p = p + 0.5 * (sum / ncount);
            update.push_back({p[0], p[1], p[2]});
        }
        else
        {
            update.push_back({p[0], p[1], p[2]});
        }
    }
    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        geometry->vertexPositions[j] = update[j] + coefficient * (geometry->vertexPositions[j] - update[j]);
    }
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}

// ======================================================================
// EXERCISE 3.2
// ========================================================================
void Exercise6::cotan_laplacian_enhance_feature()
{

    // TODO: ADD YOUR CODE HERE
    std::vector<Vector3> update;
    geometry->requireEdgeCotanWeights();
    geometry->requireVertexPositions();
    for (Vertex v : gc_mesh->vertices())
    {
        Point p = geometry->vertexPositions[v];
        if (!v.isBoundary())
        {
            double ncount = 0;
            Point sum = Vector3{0., 0., 0.};
            for (Edge e : v.adjacentEdges())
            {
                double w = geometry->edgeCotanWeights[e];
                ncount += w;
                sum += w * (geometry->vertexPositions[e.otherVertex(v)] - p);
            }
            p = p + 0.5 * (sum / ncount);
        }
        update.push_back({p[0], p[1], p[2]});
    }
    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        geometry->vertexPositions[j] = update[j] + coefficient * (geometry->vertexPositions[j] - update[j]);
    }
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}

// ======================================================================
// EXERCISE 4
// ========================================================================
void Exercise6::deform_surface(std::vector<int> &fixed_verts, const std::vector<int> &displaced_verts,
                               const Point &displacement_vector)
{

    // TODO: ADD YOUR CODE HERE
    Eigen::SparseMatrix<double> D;
    std::vector<T> DtripletList;
    Eigen::SparseMatrix<double> M;
    std::vector<T> MtripletList;
    Eigen::SparseMatrix<double> b;
    std::vector<T> btripletList;
    for (Vertex v : gc_mesh->vertices())
    {
        int i = geometry->vertexIndices[v];
        DtripletList.push_back(T(i, i, 1 / (2 * geometry->vertexDualAreas[v])));
        double totalw = 0;
        for (Edge e : v.adjacentEdges())
        {
            double w = geometry->edgeCotanWeights[e];
            totalw += geometry->edgeCotanWeights[e];
            // double w = 1;
            // totalw += 1;
            int j = geometry->vertexIndices[e.otherVertex(v)];
            MtripletList.push_back(T(i, j, w));
        }
        MtripletList.push_back(T(i, i, -1 * totalw));
    }
    // initialize b
    for (int i = 0; i < displaced_verts.size(); i++)
    {
        btripletList.push_back(T(displaced_verts[i], 0, displacement_vector[0]));
        btripletList.push_back(T(displaced_verts[i], 1, displacement_vector[1]));
        btripletList.push_back(T(displaced_verts[i], 2, displacement_vector[2]));
    }

    b.resize(DtripletList.size(), 3);
    b.setFromTriplets(btripletList.begin(), btripletList.end());
    D.resize(DtripletList.size(), DtripletList.size());
    D.setFromTriplets(DtripletList.begin(), DtripletList.end());
    M.resize(DtripletList.size(), DtripletList.size()); // should be same size as D
    M.setFromTriplets(MtripletList.begin(), MtripletList.end());

    // std::cout << "DEBUG1 " << std::endl;

    Eigen::SparseMatrix<double> L = D * M;

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    Eigen::SparseMatrix<double> Ls = L * L;

    // constrain displacement xi
    for (int i = 0; i < displaced_verts.size(); ++i)
    {
        for (int j = 0; j < Ls.cols(); ++j)
        {
            Ls.coeffRef(displaced_verts[i], j) = 0;
        }
        Ls.coeffRef(displaced_verts[i], displaced_verts[i]) = 1;
    }

    for (int i = 0; i < fixed_verts.size(); i++)
    {
        for (int j = 0; j < Ls.cols(); j++)
        {
            Ls.coeffRef(fixed_verts[i], j) = 0;
        }
        Ls.coeffRef(fixed_verts[i], fixed_verts[i]) = 1;
    }
    solver.analyzePattern(Ls);
    solver.factorize(Ls);
    Eigen::SparseMatrix<double> x = solver.solve(b);

    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        std::cout << "DEBUG2 " << Vector3{x.coeffRef(j, 0), x.coeffRef(j, 1), x.coeffRef(j, 2)} << std::endl;
        geometry->vertexPositions[j] += Vector3{x.coeffRef(j, 0), x.coeffRef(j, 1), x.coeffRef(j, 2)};
    }
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}
