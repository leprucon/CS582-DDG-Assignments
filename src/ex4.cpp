#include "include/ex4.h"
#include "include/utils.h"

// ========================================================================
// PROBLEM 1.1
// ========================================================================
void Exercise4::compute_normals_with_constant_weights()
{
    // ------------- IMPLEMENT HERE ---------
    // TODO: Compute the normals for each vertex v in the mesh using the constant
    // weights technique (see handout) and store inside vertex_normals
    // ------------- IMPLEMENT HERE ---------
    vertex_normals.clear(); //???
    geometry->requireFaceNormals();
    for (Vertex v : gc_mesh->vertices())
    {
        Point normalSum{0., 0., 0.};
        for (Face f : v.adjacentFaces())
        {
            double x = (geometry->faceNormals[f].normalize()[0]); // not sure if the face normals in geometry central are normalized by default, so I did it just to be safe.
            double y = (geometry->faceNormals[f].normalize()[1]);
            double z = (geometry->faceNormals[f].normalize()[2]);
            normalSum = Vector3{normalSum[0] + x, normalSum[1] + y, normalSum[2] + z};
        }
        normalSum = normalSum.normalize();
        vertex_normals.push_back({normalSum[0], normalSum[1], normalSum[2]});
    }
}

// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Exercise4::compute_normals_by_area_weights()
{
    // ------------- IMPLEMENT HERE ---------
    // TODO: Compute the normals for each vertex v in the mesh using weights
    // proportional to the dual area (see handout)
    // and store inside vertex_normals
    // ------------- IMPLEMENT HERE ---------
    vertex_normals.clear(); //???
    geometry->requireFaceNormals();
    geometry->requireFaceAreas();
    for (Vertex v : gc_mesh->vertices())
    {
        Point normalSum{0., 0., 0.};
        for (Face f : v.adjacentFaces())
        {
            double x = (geometry->faceAreas[f]) * (geometry->faceNormals[f].normalize()[0]); // not sure if the face normals in geometry central are normalized by default, so I did it just to be safe.
            double y = (geometry->faceAreas[f]) * (geometry->faceNormals[f].normalize()[1]);
            double z = (geometry->faceAreas[f]) * (geometry->faceNormals[f].normalize()[2]);
            normalSum = Vector3{normalSum[0] + x, normalSum[1] + y, normalSum[2] + z};
        }
        normalSum = normalSum.normalize();
        vertex_normals.push_back({normalSum[0], normalSum[1], normalSum[2]});
    }
}

// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Exercise4::compute_normals_with_angle_weights()
{
    // ------------- IMPLEMENT HERE ---------
    // TODO: Compute the normals for each vertex v in the mesh using weights
    // proportional to incident angles (see handout)
    // and store inside vertex_normals
    // ------------- IMPLEMENT HERE ---------
    // vertex_normals.clear(); //???
    // geometry->requireFaceNormals();
    // geometry->requireFaceAreas();
    // for (Vertex v : gc_mesh->vertices())
    // {
    //     Point normalSum{0., 0., 0.};
    //     for (Face f : v.adjacentFaces())
    //     {
    //         double x = (geometry->faceAreas[f]) * (geometry->faceNormals[f].normalize()[0]); // not sure if the face normals in geometry central are normalized by default, so I did it just to be safe.
    //         double y = (geometry->faceAreas[f]) * (geometry->faceNormals[f].normalize()[1]);
    //         double z = (geometry->faceAreas[f]) * (geometry->faceNormals[f].normalize()[2]);
    //         normalSum = Vector3{normalSum[0] + x, normalSum[1] + y, normalSum[2] + z};
    //     }
    //     normalSum = normalSum.normalize();
    //     vertex_normals.push_back({normalSum[0], normalSum[1], normalSum[2]});
    // }
    vertex_normals.clear();
    geometry->requireVertexNormals();
    for (Vertex v : gc_mesh->vertices())
    {
        Point normalSum = geometry->vertexNormals[v];
        normalSum = normalSum.normalize();
        vertex_normals.push_back({normalSum[0], normalSum[1], normalSum[2]});
    }
}

// ========================================================================
// EXERCISE 2.1
// ========================================================================
void Exercise4::calc_uniform_laplacian()
{
    // // ------------- IMPLEMENT HERE ---------
    // // TODO:
    // // For each non-boundary vertex, compute uniform Laplacian operator vector
    // // and store the vector length in the vertex property of the
    // // mesh called vertex_curvature_weights.
    // // Store min and max values in min_curvature and max_curvature.
    // // ------------- IMPLEMENT HERE ---------
    vertex_curvature_weights.clear();
    geometry->requireVertexPositions();
    max_curvature = -DBL_MAX;
    min_curvature = DBL_MAX;
    int vertexindex = 0;
    for (Vertex v : gc_mesh->vertices())
    {
        double x = 0.;
        double y = 0.;
        double z = 0.;
        int n = 0;
        for (Vertex vi : v.adjacentVertices())
        {
            x += geometry->vertexPositions[vi][0] - geometry->vertexPositions[v][0];
            y += geometry->vertexPositions[vi][1] - geometry->vertexPositions[v][1];
            z += geometry->vertexPositions[vi][2] - geometry->vertexPositions[v][2];
            n++;
        }
        if (v.isBoundary())
        {
            vertex_curvature_weights.push_back(0); // placeholder weight?
        }
        else
        {
            x = x / n;
            y = y / n;
            z = z / n;
            // double len = std::pow((std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)), 0.5);                                           // not sure what 'n' is in the homework description... just used vector length like described in this file's TODO.
            double len = (vertex_normals[vertexindex][0] * x + vertex_normals[vertexindex][1] * y + vertex_normals[vertexindex][2] * z) / 2; // assuming n is the normal of the vector?
            vertex_curvature_weights.push_back(len);
            if (len > max_curvature)
            {
                max_curvature = len;
            }
            if (len < min_curvature)
            {
                min_curvature = len;
            }
        }
        vertexindex++;
    }

    std::cout << "Min. uniform value is: " << min_curvature << std::endl;
    std::cout << "Max. uniform value is: " << max_curvature << std::endl;
}

// ========================================================================
// EXERCISE 2.2
// ========================================================================
void Exercise4::calc_mean_curvature()
{
    // ------------- IMPLEMENT HERE ---------
    // // TODO:
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // ------------- IMPLEMENT HERE ---------
    vertex_curvature_weights.clear();
    geometry->requireVertexPositions();
    max_curvature = -DBL_MAX;
    min_curvature = DBL_MAX;
    int vertexindex = 0;
    for (Vertex v : gc_mesh->vertices())
    {
        double x = 0.;
        double y = 0.;
        double z = 0.;
        for (Edge e : v.adjacentEdges())
        {
            Vertex vi;
            if (e.adjacentVertices()[0] == v)
            {
                vi = e.adjacentVertices()[1];
            }
            else
            {
                vi = e.adjacentVertices()[0];
            }
            x += geometry->edgeCotanWeights[e] * (geometry->vertexPositions[vi][0] - geometry->vertexPositions[v][0]);
            y += geometry->edgeCotanWeights[e] * (geometry->vertexPositions[vi][1] - geometry->vertexPositions[v][1]);
            z += geometry->edgeCotanWeights[e] * (geometry->vertexPositions[vi][2] - geometry->vertexPositions[v][2]);
        }
        if (v.isBoundary())
        {
            vertex_curvature_weights.push_back(0); // placeholder weight?
        }
        else
        {
            double A = geometry->vertexDualAreas[v];
            x = x / (2 * A);
            y = y / (2 * A);
            z = z / (2 * A);
            // double len = std::pow((std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)), 0.5);                                           // not sure what 'n' is in the homework description... just used vector length like described in this file's TODO.
            double len = (vertex_normals[vertexindex][0] * x + vertex_normals[vertexindex][1] * y + vertex_normals[vertexindex][2] * z) / 2; // assuming n is the normal of the vector?
            vertex_curvature_weights.push_back(len);
            if (len > max_curvature)
            {
                max_curvature = len;
            }
            if (len < min_curvature)
            {
                min_curvature = len;
            }
        }
        vertexindex++;
    }

    std::cout << "Min. Laplace-Beltrami value is: " << min_curvature << std::endl;
    std::cout << "Max. Laplace-Beltrami value is: " << max_curvature << std::endl;
}

// ========================================================================
// EXERCISE 2.3
// ========================================================================
void Exercise4::calc_gauss_curvature()
{
    // ------------- IMPLEMENT HERE ---------
    // // TODO:
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of dot products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the vertex_weight_ property for the area weight.

    // ------------- IMPLEMENT HERE ---------
    // IM ASSUMING THAT THE TODO TEXT IS OUTDATED
    vertex_curvature_weights.clear();
    max_curvature = -DBL_MAX;
    min_curvature = DBL_MAX;
    double twopi = 3.14159265358979323846 * 2;
    for (Vertex v : gc_mesh->vertices())
    {
        if (v.isBoundary())
        {
            vertex_curvature_weights.push_back(0); // placeholder weight?
        }
        else
        {
            double A = geometry->vertexDualAreas[v];
            double anglesum = geometry->vertexAngleSums[v];
            double G = (twopi - anglesum) / A;
            // double len = std::pow((std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2)), 0.5);                                           // not sure what 'n' is in the homework description... just used vector length like described in this file's TODO.
            vertex_curvature_weights.push_back(G);
            if (G > max_curvature)
            {
                max_curvature = G;
            }
            if (G < min_curvature)
            {
                min_curvature = G;
            }
        }
    }
    std::cout << "Min. Gaussian value is: " << min_curvature << std::endl;
    std::cout << "Max. Gaussian value is: " << max_curvature << std::endl;
}

//====================================================================================================================//

void Exercise4::show_normal()
{
    switch (normal_type)
    {
    case 0:
        compute_normals_with_constant_weights();
        break;
    case 1:
        compute_normals_by_area_weights();
        break;
    case 2:
        compute_normals_with_angle_weights();
        break;
    }
    auto vvq = ps_mesh->addVertexVectorQuantity("normal", vertex_normals);
    vvq->setEnabled(true);
}

void Exercise4::show_curvature()
{
    switch (curvature_type)
    {
    case 0:
        calc_uniform_laplacian();
        break;
    case 1:
        calc_mean_curvature();
        break;
    case 2:
        calc_gauss_curvature();
        break;
    }
    color_coding_d(vertex_curvature_weights);
}

void Exercise4::color_coding_d(const std::vector<double> data, const double _min_value,
                               const double _max_value, const int _bound)
{
    std::vector<double> datavalue(data.size());
    for (size_t i = 0; i < data.size(); i++)
    {
        datavalue[i] = data[i];
    }
    auto min_value = _min_value;
    auto max_value = _max_value;

    if (min_value == 0 && max_value == 0)
    {
        // discard upper and lower bound
        auto n = data.size() - 1;
        auto i = n / _bound;

        std::sort(datavalue.begin(), datavalue.end());
        min_value = datavalue[i];
        max_value = datavalue[n - 1 - i];
    }

    const auto range = max_value - min_value;
    std::vector<std::array<double, 3>> psColor;

    for (auto vh : gc_mesh->vertices())
    {
        auto t = ((data[vh.getIndex()]) - min_value) / range;
        auto color = util::map_val2color(t, 0, 1.0);
        psColor.push_back({color.x, color.y, color.z});
    }
    auto vcq = ps_mesh->addVertexColorQuantity("Vertex Colors", psColor);
    vcq->setEnabled(true);
}
