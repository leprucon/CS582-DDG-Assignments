#include "include/ex7.h"
#include <iostream>
#include "geometrycentral/surface/halfedge_element_types.h"

// ======================================================================
// EXERCISE 1.1 Mapping to circle
// ========================================================================
void Exercise7::map_suface_boundary_to_circle()
{

    // initialize all points to origin

    for (Vertex v : tex_mesh->vertices())
    {
        tex_geometry->vertexPositions[v] = _origin;
    }

    // TODO: ADD YOUR CODE HERE
    int validcount = 0;
    BoundaryLoop b;
    for (BoundaryLoop bi : tex_mesh->boundaryLoops())
    {
        b = bi;
        validcount++;
    }
    if (validcount > 1)
    {
        std::cout << "Not a disk topology!" << std::endl;
    }
    // double _circumference = 2 * _radius * 3.141592653589;
    double total_boundary_length = 0;
    int numBoundVert = tex_mesh->nVertices() - tex_mesh->nInteriorVertices();
    for (Vertex v : b.adjacentVertices())
    {
        Vector3 diff = geometry->vertexPositions[v.getIndex()] - geometry->vertexPositions[(v.getIndex() + 1) % (numBoundVert)];
        total_boundary_length += std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.y * diff.y);
    }
    double progress = 0;
    for (Vertex v : b.adjacentVertices())
    {
        Vector3 diff = geometry->vertexPositions[v.getIndex()] - geometry->vertexPositions[(v.getIndex() + 1) % (numBoundVert)];
        double len = std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.y * diff.y);
        progress = progress + ((len / total_boundary_length) * 2 * 3.14159265358979);
        tex_geometry->vertexPositions[v] = _origin + Vector3{_radius * cos(progress), _radius * sin(progress), 0};
    }

    // // assuming that the selected mesh has a disk topology... (does not check)
    // Vertex current;
    // Vertex start;
    // for (Vertex v : gc_mesh->vertices())
    // {
    //     if (v.isBoundary())
    //     {
    //         current = v;
    //         start = current;
    //         break;
    //     }
    // }
    // Vertex prev = current;
    // int stop = 1;
    // std::vector<Vertex> boundverts;
    // // geometry->requireVertexIndices();
    // while (true)
    // {
    //     if (current == start)
    //     {
    //         if (stop == 0)
    //         {
    //             break;
    //         }
    //         stop--;
    //     }
    //     boundverts.push_back(current);
    //     for (Edge en : current.adjacentEdges())
    //     {
    //         Vertex vn = en.otherVertex(current);
    //         if (en.isBoundary() && vn != prev)
    //         {
    //             prev = current;
    //             current = vn;
    //             break;
    //         }
    //     }
    // }

    // double totallen = 0;
    // bool first = true;
    // Vertex back = boundverts.back();
    // for (Vertex v : boundverts) // compute total boundary length
    // {
    //     if (first)
    //     {
    //         first = false;
    //         Vector3 diff = geometry->vertexPositions[v] - geometry->vertexPositions[back];
    //         totallen += std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.y * diff.y);
    //         back = v;
    //     }
    //     else
    //     {
    //         Vector3 diff = geometry->vertexPositions[v] - geometry->vertexPositions[back];
    //         totallen += std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.y * diff.y);
    //     }
    // }

    // first = true;
    // back = boundverts.back();
    // double progress = 0;
    // for (Vertex v : boundverts) // distribute boundary vectors along circle
    // {
    //     if (first)
    //     {
    //         first = false;
    //         tex_geometry->vertexPositions[v] = _origin + Vector3{_radius * cos(progress), _radius * sin(progress), 0};
    //         back = v;
    //     }
    //     else
    //     {
    //         Vector3 diff = geometry->vertexPositions[v] - geometry->vertexPositions[back];
    //         double len = std::sqrt(diff.x * diff.x + diff.y * diff.y + diff.y * diff.y);
    //         progress = progress + ((len / totallen) * 2 * 3.14159265358979);
    //         tex_geometry->vertexPositions[v] = _origin + Vector3{_radius * cos(progress), _radius * sin(progress), 0};
    //     }
    // }
    // std::cout << "PROGRESS: " << (progress / 2) << std::endl;
    // // tex_geometry->vertexPositions[v] = _origin + Vector3{cos(0), sin(0), 0};

    // // for (Vertex v : gc_mesh->vertices()) // distribute boundary vectors along circle
    // // {
    // //     if (v.isBoundary())
    // //     {
    // //         tex_geometry->vertexPositions[v] = geometry->vertexPositions[v];
    // //     }
    // // }

    // register mesh with polyscope
    tex_geometry->refreshQuantities();
    ps_tex_mesh = polyscope::registerSurfaceMesh("tex_coords_circle", tex_geometry->vertexPositions, tex_mesh->getFaceVertexList());
    ps_tex_mesh->setEdgeWidth(1.);
    ps_tex_mesh->setSurfaceColor({1., 1., 1.});
    is_mapped = true;
}

// ======================================================================
// EXERCISE 1.2 Interactively smoothing the texture
// ========================================================================
void Exercise7::explicitly_smooth_texture()
{
    // TODO: ADD YOUR CODE HERE
    for (int i = 0; i < num_iter; ++i)
    {
        geometry->requireEdgeCotanWeights();
        for (Vertex v : tex_mesh->vertices())
        {
            if (!v.isBoundary())
            {
                double weightsum = 0;
                Vector3 possum = {0., 0., 0.};
                geometry->requireHalfedgeCotanWeights();
                geometry->requireCornerAngles();
                for (Edge e : v.adjacentEdges())
                {
                    double rectifiedWeight = geometry->edgeCotanWeight(e); // I am the dumbest mf in the history of man. I spent two days stuck on this problem because I used tex_geometry here instead of geometry.
                    if (rectifiedWeight <= 0)
                    {
                        // Halfedge he = e.halfedge();
                        // double a = geometry->cornerAngle(he.next().corner()) + geometry->cornerAngle(he.twin().corner());
                        // double b = geometry->cornerAngle(he.corner()) + geometry->cornerAngle(he.twin().next().corner());
                        // rectifiedWeight = (0.5 / std::tan(a)) + (0.5 / std::tan(b)); // pretend it's delaunay triangulation by edge flipping???

                        rectifiedWeight = 0;
                    }
                    if (rectifiedWeight != rectifiedWeight || std::isinf(rectifiedWeight))
                    {
                        // possum += (tex_geometry->vertexPositions[e.otherVertex(v)]);
                        // weightsum += 1;
                        // continue;
                        rectifiedWeight = 0; // DBL_MAX?
                    }
                    // std::cout << "Weight: " << rectifiedWeight << std::endl;
                    possum += (tex_geometry->vertexPositions[e.otherVertex(v)]) * (rectifiedWeight);
                    weightsum += rectifiedWeight;
                    // std::cout << "Weight: " << tex_geometry->edgeCotanWeight(e) << std::endl;
                    // possum += (tex_geometry->vertexPositions[e.otherVertex(v)]);
                    // weightsum += 1;
                }
                if (weightsum == 0)
                {
                    std::cout << "ZERO" << std::endl;
                }
                else
                {
                    tex_geometry->vertexPositions[v] = possum / weightsum;
                    _tex_coords[v.getIndex()] = {tex_geometry->vertexPositions[v].x, tex_geometry->vertexPositions[v].y};
                }
                // std::cout << "SUM: " << weightsum << std::endl;
            }
        }

        geometry->refreshQuantities();
        auto qParam = ps_mesh->addVertexParameterizationQuantity("tex_coords", _tex_coords);
        qParam->setCheckerSize(1.);
        qParam->setEnabled(true);
        ps_tex_mesh->updateVertexPositions(tex_geometry->vertexPositions);
    }
}

// ======================================================================
// EXERCISE 1.3 Implicitly smoothing the texture
// ========================================================================
void Exercise7::implicitly_smooth_texture()
{

    // TODO: ADD YOUR CODE HERE
    geometry->requireVertexIndices();
    tex_geometry->requireVertexIndices();
    Vertex geometry_getVertex[gc_mesh->nVertices()];
    Vertex tex_geometry_getVertex[tex_mesh->nVertices()];
    int sortedindex = 0;
    for (Vertex v : gc_mesh->vertices())
    {
        if (!v.isBoundary())
        {
            geometry->vertexIndices[v] = sortedindex;
            geometry_getVertex[sortedindex] = v;
            sortedindex++;
        }
    }
    for (Vertex v : gc_mesh->vertices())
    {
        if (v.isBoundary())
        {
            geometry->vertexIndices[v] = sortedindex;
            geometry_getVertex[sortedindex] = v;
            sortedindex++;
        }
    }
    sortedindex = 0;
    for (Vertex v : tex_mesh->vertices()) // Should be unnecessary, but I was running into a bunch of bugs that I couldn't remotely understand so I'm just adding code to be as thorough as possible.
    {
        if (!v.isBoundary())
        {
            tex_geometry->vertexIndices[v] = sortedindex;
            tex_geometry_getVertex[sortedindex] = v;
            sortedindex++;
        }
    }
    for (Vertex v : tex_mesh->vertices())
    {
        if (v.isBoundary())
        {
            tex_geometry->vertexIndices[v] = sortedindex;
            tex_geometry_getVertex[sortedindex] = v;
            sortedindex++;
        }
    }
    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> bneg; // negative b
    std::vector<T> MtripletList;
    std::vector<T> bnegtripletList;
    int sizen = 0; // could just use geometry central quantity?
    for (Vertex v : gc_mesh->vertices())
    {
        if (!v.isBoundary())
        {
            sizen++;
            int i = geometry->vertexIndices[v];
            double totalw = 0;
            Vector3 cotansum = Vector3{0., 0., 0.};
            for (Edge e : v.adjacentEdges())
            {
                double w = geometry->edgeCotanWeights[e];
                if (!e.otherVertex(v).isBoundary())
                {
                    int j = geometry->vertexIndices[e.otherVertex(v)];
                    MtripletList.push_back(T(i, j, w));
                    totalw += geometry->edgeCotanWeights[e];
                }
                else
                {
                    totalw += geometry->edgeCotanWeights[e];
                    cotansum += w * tex_geometry->vertexPositions[e.otherVertex(v).getIndex()];
                }
            }
            MtripletList.push_back(T(i, i, -1 * totalw));
            bnegtripletList.push_back(T(i, 0, cotansum.x));
            bnegtripletList.push_back(T(i, 1, cotansum.y));
            bnegtripletList.push_back(T(i, 2, cotansum.z));
        }
    }
    M.resize(sizen, sizen);
    M.setFromTriplets(MtripletList.begin(), MtripletList.end());
    bneg.resize(sizen, 3);
    bneg.setFromTriplets(bnegtripletList.begin(), bnegtripletList.end());

    // geometry->requireCotanLaplacian(); //TEST
    // Eigen::LLT<Eigen::MatrixXd> A_llt(-1 * M);
    // if (A_llt.info() == Eigen::NumericalIssue)
    // {
    //     std::cout << "BAD" << std::endl;
    // }
    // std::cout << "M: " << -1 * M << std::endl;

    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    // solver.analyzePattern(-1 * M);
    // solver.factorize(-1 * M);
    // solver.compute(-1 * M);
    // Eigen::MatrixXd update = solver.solve(bneg);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(-1 * M);
    solver.factorize(-1 * M);
    solver.compute(-1 * M);
    Eigen::MatrixXd update = solver.solve(bneg);

    // std::cout << "RESULT: " << update << std::endl;

    int skip = 0;
    for (size_t j = 0; j < tex_mesh->nVertices(); j++)
    {
        if (tex_geometry_getVertex[j].isBoundary())
        {
            skip++;
            continue;
        }
        int jskip = j - skip;
        tex_geometry->vertexPositions[tex_geometry_getVertex[jskip]] = {update(jskip, 0), update(jskip, 1), update(jskip, 2)};
        _tex_coords[tex_geometry_getVertex[jskip].getIndex()] = {update(jskip, 0), update(jskip, 1)};
    }

    geometry->refreshQuantities();
    auto qParam = ps_mesh->addVertexParameterizationQuantity("tex_coords", _tex_coords);
    qParam->setCheckerSize(1.);
    qParam->setEnabled(true);
    ps_tex_mesh->updateVertexPositions(tex_geometry->vertexPositions);
}

// ======================================================================
// EXERCISE 2 Minimal Surfaces
// ======================================================================
void Exercise7::minimal_surface()
{

    // TODO: ADD YOUR CODE HERE
    Eigen::SparseMatrix<double> M;
    std::vector<T> MtripletList;
    Eigen::SparseMatrix<double> Zero;
    std::vector<T> ZerotripletList;
    for (Vertex v : gc_mesh->vertices())
    {
        int i = v.getIndex();
        if (!v.isBoundary())
        {
            double totalw = 0;
            for (Edge e : v.adjacentEdges())
            {
                double w = geometry->edgeCotanWeights[e];
                totalw += geometry->edgeCotanWeights[e];
                // double w = 1;
                // totalw += 1;
                int j = e.otherVertex(v).getIndex();
                MtripletList.push_back(T(i, j, w));
            }
            MtripletList.push_back(T(i, i, -1 * totalw));
        }
        else
        {
            MtripletList.push_back(T(i, i, 1));
            ZerotripletList.push_back(T(v.getIndex(), 0, geometry->vertexPositions[v].x));
            ZerotripletList.push_back(T(v.getIndex(), 1, geometry->vertexPositions[v].y));
            ZerotripletList.push_back(T(v.getIndex(), 2, geometry->vertexPositions[v].z));
        }
    }
    M.resize(gc_mesh->nVertices(), gc_mesh->nVertices()); // should be same size as D
    M.setFromTriplets(MtripletList.begin(), MtripletList.end());
    Zero.resize(gc_mesh->nVertices(), 3);
    Zero.setFromTriplets(ZerotripletList.begin(), ZerotripletList.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::SparseMatrix<double> x = solver.solve(Zero);

    for (size_t j = 0; j < gc_mesh->nVertices(); j++)
    {
        geometry->vertexPositions[j] = Vector3{x.coeffRef(j, 0), x.coeffRef(j, 1), x.coeffRef(j, 2)};
    }

    geometry->refreshQuantities();
    ps_mesh->updateVertexPositions(geometry->vertexPositions);
}