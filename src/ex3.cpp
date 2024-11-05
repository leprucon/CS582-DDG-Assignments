#include "include/ex3.h"

void Exercise3::generate_curve_simple()
{
    std::cout << "Generating simple curve" << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    Point center = {0., 0., 0.};
    double radius = 1;

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    first = true;
    for (size_t i = 0; i < num_vertices; ++i)
    {
        Point pt;
        double frac = static_cast<double>(i) / static_cast<double>(num_vertices);
        pt[0] = center[0] + radius * (cos(2. * M_PI * frac) + distribution(generator));
        pt[1] = center[0] + radius * (sin(2. * M_PI * frac) + distribution(generator));
        pt[2] = 0.;
        points.push_back(pt);
    }

    for (size_t i = 0; i < num_vertices; ++i)
    {
        edges.push_back({i, (i + 1) % num_vertices});
    }
    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::generate_curve_figure_eight()
{
    std::cout << "Generating figure-8" << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    // TODO: ADD YOUR CODE HERE
    // y^2 - x^2(a^2 - x^2)
    // y^2 = sin(x^2)
    first = true;
    for (size_t i = 0; i <= num_vertices / 2; ++i)
    {
        Point pt;
        double frac = static_cast<double>(i) / (static_cast<double>(num_vertices) / 2);
        double x = frac * 2 - 1;
        pt[0] = x + distribution(generator);
        if (x <= 0)
        {
            pt[1] = pow(pow(x, 2) - pow(x, 4), 0.5) + distribution(generator);
        }
        if (x > 0)
        {
            pt[1] = pow(pow(x, 2) - pow(x, 4), 0.5) * -1 + distribution(generator);
        }
        pt[2] = 0.;
        points.push_back(pt);
    }
    for (size_t i = num_vertices - 1; i > num_vertices / 2; --i)
    {
        Point pt;
        double frac = static_cast<double>(i - (num_vertices / 2)) / (static_cast<double>(num_vertices) / 2);
        double x = frac * 2 - 1;
        pt[0] = x + distribution(generator);
        if (x <= 0)
        {
            pt[1] = pow(pow(x, 2) - pow(x, 4), 0.5) * -1 + distribution(generator);
        }
        if (x > 0)
        {
            pt[1] = pow(pow(x, 2) - pow(x, 4), 0.5) + distribution(generator);
        }
        pt[2] = 0.;
        points.push_back(pt);
    }

    for (size_t i = 0; i < num_vertices; ++i)
    {
        edges.push_back({i, (i + 1) % num_vertices});
    }
    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::generate_curve_limacon()
{
    std::cout << "Generating limaçon" << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    // TODO: ADD YOUR CODE HERE
    first = true;
    for (size_t i = 0; i < num_vertices; ++i)
    {
        Point pt;
        double frac = static_cast<double>(i) / static_cast<double>(num_vertices);
        double radius = 0.5 + 2 * cos(2. * M_PI * frac);
        pt[0] = radius * cos(2. * M_PI * frac) + distribution(generator);
        pt[1] = radius * sin(2. * M_PI * frac) + distribution(generator);
        pt[2] = 0.;
        points.push_back(pt);
    }

    for (size_t i = 0; i < num_vertices; ++i)
    {
        edges.push_back({i, (i + 1) % num_vertices});
    }
    Curve = polyscope::registerCurveNetwork("curve", points, edges);
}

void Exercise3::generate_curve_3d()
{
    std::cout << "Generating 3D curve" << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    points.clear();
    points.reserve(num_vertices);
    edges.clear();
    edges.reserve(num_vertices);

    // TODO: ADD YOUR CODE HERE
    first = true;
    double radius = 1;
    for (size_t i = 0; i < num_vertices; ++i)
    {
        Point pt;
        double frac = static_cast<double>(i) / static_cast<double>(num_vertices);
        pt[0] = radius * (cos(2. * M_PI * frac) + distribution(generator));
        pt[1] = radius * (sin(2. * M_PI * frac) + distribution(generator));
        pt[2] = -1 * radius * (sin(2. * M_PI * frac) + distribution(generator));
        points.push_back(pt);
    }

    for (size_t i = 0; i < num_vertices; ++i)
    {
        edges.push_back({i, (i + 1) % num_vertices});
    }
    Curve = polyscope::registerCurveNetwork("curve", points, edges);
    // generate_world_axes(); // uncomment to visualize world axes
}

void Exercise3::laplacian_smoothing()
{

    // TODO: ADD YOUR CODE HERE
    for (size_t q = 0; q < num_iter; ++q)
    {
        if (first)
        {
            maxX = points[0][0];
            minX = points[0][0];
            maxY = points[0][1];
            minY = points[0][1];
            maxZ = points[0][2];
            minZ = points[0][2];
        }
        for (size_t i = 0; i < num_vertices; i++)
        {
            if (first)
            {
                if (maxX < points[i][0])
                {
                    maxX = points[i][0];
                }
                if (minX > points[i][0])
                {
                    minX = points[i][0];
                }
                if (maxY < points[i][1])
                {
                    maxY = points[i][1];
                }
                if (minY > points[i][1])
                {
                    minY = points[i][1];
                }
                if (maxZ < points[i][2])
                {
                    maxZ = points[i][2];
                }
                if (minZ > points[i][2])
                {
                    minZ = points[i][2];
                }
            }
            int wrap = static_cast<int>(num_vertices);
            points[(i + 1) % wrap][0] = (1 - epsilon) * points[(i + 1) % wrap][0] + epsilon * (points[(i) % wrap][0] + points[(i + 2) % wrap][0]) / 2;
            points[(i + 1) % wrap][1] = (1 - epsilon) * points[(i + 1) % wrap][1] + epsilon * (points[(i) % wrap][1] + points[(i + 2) % wrap][1]) / 2;
            points[(i + 1) % wrap][2] = (1 - epsilon) * points[(i + 1) % wrap][2] + epsilon * (points[(i) % wrap][2] + points[(i + 2) % wrap][2]) / 2;
        }
        first = false;
        double maxX2 = points[0][0];
        double minX2 = points[0][0];
        double maxY2 = points[0][1];
        double minY2 = points[0][1];
        double maxZ2 = points[0][2];
        double minZ2 = points[0][2];
        double centx;
        double centy;
        double centz;
        for (size_t i = 0; i < num_vertices; i++)
        {
            if (maxX2 < points[i][0])
            {
                maxX2 = points[i][0];
            }
            if (minX2 > points[i][0])
            {
                minX2 = points[i][0];
            }
            if (maxY2 < points[i][1])
            {
                maxY2 = points[i][1];
            }
            if (minY2 > points[i][1])
            {
                minY2 = points[i][1];
            }
            if (maxZ2 < points[i][2])
            {
                maxZ2 = points[i][2];
            }
            if (minZ2 > points[i][2])
            {
                minZ2 = points[i][2];
            }
            centx += points[i][0];
            centy += points[i][1];
            centz += points[i][2];
        }
        centx = centx / num_vertices;
        centy = centy / num_vertices;
        centz = centz / num_vertices;

        for (size_t i = 0; i < num_vertices; i++)
        {
            points[i][0] = (points[i][0] - centx) * ((maxX - minX) / (maxX2 - minX2)) + centx;
            points[i][1] = (points[i][1] - centy) * ((maxY - minY) / (maxY2 - minY2)) + centy;
            // points[i][2] = (points[i][2] - centz) * ((maxZ - minZ) / (maxZ2 - minZ2)) + centz; // divide by zero edge case
        }
    }
    Curve->updateNodePositions(points);
}

void Exercise3::osculating_circle()
{

    // TODO: ADD YOUR CODE HERE
    for (size_t q = 0; q < num_iter; ++q)
    {
        if (first)
        {
            maxX = points[0][0];
            minX = points[0][0];
            maxY = points[0][1];
            minY = points[0][1];
            maxZ = points[0][2];
            minZ = points[0][2];
        }
        for (size_t i = 0; i < num_vertices; i++)
        {
            if (first)
            {
                if (maxX < points[i][0])
                {
                    maxX = points[i][0];
                }
                if (minX > points[i][0])
                {
                    minX = points[i][0];
                }
                if (maxY < points[i][1])
                {
                    maxY = points[i][1];
                }
                if (minY > points[i][1])
                {
                    minY = points[i][1];
                }
                if (maxZ < points[i][2])
                {
                    maxZ = points[i][2];
                }
                if (minZ > points[i][2])
                {
                    minZ = points[i][2];
                }
            }
            int wrap = static_cast<int>(num_vertices);
            Point circCircCent = CircumscribedCircle(points[(i) % wrap], points[(i + 1) % wrap], points[(i + 2) % wrap]);
            glm::vec3 mag = {circCircCent[0] - points[(i + 1) % wrap][0],
                             circCircCent[1] - points[(i + 1) % wrap][1],
                             circCircCent[2] - points[(i + 1) % wrap][2]};
            double len2 = pow(glm::length(mag), 2);
            points[(i + 1) % wrap][0] = points[(i + 1) % wrap][0] + epsilon * (circCircCent[0] - points[(i + 1) % wrap][0]) / (len2);
            points[(i + 1) % wrap][1] = points[(i + 1) % wrap][1] + epsilon * (circCircCent[1] - points[(i + 1) % wrap][1]) / (len2);
            points[(i + 1) % wrap][2] = points[(i + 1) % wrap][2] + epsilon * (circCircCent[2] - points[(i + 1) % wrap][2]) / (len2);
        }
        first = false;
        double maxX2 = points[0][0];
        double minX2 = points[0][0];
        double maxY2 = points[0][1];
        double minY2 = points[0][1];
        double maxZ2 = points[0][2];
        double minZ2 = points[0][2];
        double centx;
        double centy;
        double centz;
        for (size_t i = 0; i < num_vertices; i++)
        {
            if (maxX2 < points[i][0])
            {
                maxX2 = points[i][0];
            }
            if (minX2 > points[i][0])
            {
                minX2 = points[i][0];
            }
            if (maxY2 < points[i][1])
            {
                maxY2 = points[i][1];
            }
            if (minY2 > points[i][1])
            {
                minY2 = points[i][1];
            }
            if (maxZ2 < points[i][2])
            {
                maxZ2 = points[i][2];
            }
            if (minZ2 > points[i][2])
            {
                minZ2 = points[i][2];
            }
            centx += points[i][0];
            centy += points[i][1];
            centz += points[i][2];
        }
        centx = centx / num_vertices;
        centy = centy / num_vertices;
        centz = centz / num_vertices;

        for (size_t i = 0; i < num_vertices; i++)
        {
            points[i][0] = (points[i][0] - centx) * ((maxX - minX) / (maxX2 - minX2)) + centx;
            points[i][1] = (points[i][1] - centy) * ((maxY - minY) / (maxY2 - minY2)) + centy;
            // points[i][2] = (points[i][2] - centz) * ((maxZ - minZ) / (maxZ2 - minZ2)) + centz; // divide by zero edge case
        }
    }
    Curve->updateNodePositions(points);
}

glm::vec3 Exercise3::CircumscribedCircle(const Point &a, const Point &b, const Point &c)
{
    Point center = {0., 0., 0.};

    // TODO: ADD YOUR CODE HERE
    // calculate circumcenter of triangle?
    // taken from gamedev.stackexchange.com/questions/60630
    // glm::vec3 ac = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    // glm::vec3 ab = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    // glm::vec3 abxac = glm::cross(ab, ac);
    // glm::vec3 abxacXab = glm::cross(abxac, ab);
    // glm::vec3 acXabxac = glm::cross(ac, abxac);
    // double cx = a[0] + (abxacXab[0] * pow(glm::length(ac), 2) + acXabxac[0] * pow(glm::length(ac), 2)) / (2 * pow(glm::length(abxac), 2));
    // double cy = a[1] + (abxacXab[1] * pow(glm::length(ac), 2) + acXabxac[1] * pow(glm::length(ac), 2)) / (2 * pow(glm::length(abxac), 2));
    // double cz = a[2] + (abxacXab[2] * pow(glm::length(ac), 2) + acXabxac[2] * pow(glm::length(ac), 2)) / (2 * pow(glm::length(abxac), 2));
    // center = {cx, cy, cz};
    double a2[3] = {a[0], a[1], a[2]};
    double b2[3] = {b[0], b[1], b[2]};
    double c2[3] = {c[0], c[1], c[2]};
    double center2[3];
    tricircumcenter3d(a2, b2, c2, center2);
    center[0] = center2[0] + a2[0];
    center[1] += center2[1] + a2[1];
    center[2] += center2[2] + a2[2];
    std::cout << "INPUTS: " << a << ", " << b << ", " << c << std::endl;
    std::cout << "RESULT: " << center << std::endl;
    return center;
}

void Exercise3::generate_world_axes()
{
    PCPointCloud *ax_x = new PCPointCloud(1);
    PointPositionGeometry *geo_x = new PointPositionGeometry(*ax_x);
    geo_x->positions[0] = {0., 0., -1.};
    PointData<Point> ax_x_q(*ax_x);
    ax_x_q[0] = {0.25, 0., 0.};
    auto pc_x = polyscope::registerPointCloud("World x-axis", geo_x->positions);
    pc_x->setPointColor({0., 0., 0.});
    auto vis_x = pc_x->addVectorQuantity("Axes", ax_x_q, VectorType::AMBIENT);
    vis_x->setVectorColor({1., 0., 0.});
    vis_x->setVectorRadius(0.01);
    vis_x->setEnabled(true);

    PCPointCloud *ax_y = new PCPointCloud(1);
    PointPositionGeometry *geo_y = new PointPositionGeometry(*ax_y);
    geo_y->positions[0] = {0., 0., -1.};
    PointData<Point> ax_y_q(*ax_y);
    ax_y_q[0] = {0., 0.25, 0.};
    auto pc_y = polyscope::registerPointCloud("World y-axis", geo_y->positions);
    pc_y->setPointColor({0., 0., 0.});
    auto vis_y = pc_y->addVectorQuantity("Axes", ax_y_q, VectorType::AMBIENT);
    vis_y->setVectorColor({0., 1., 0.});
    vis_y->setVectorRadius(0.01);
    vis_y->setEnabled(true);

    PCPointCloud *ax_z = new PCPointCloud(1);
    PointPositionGeometry *geo_z = new PointPositionGeometry(*ax_z);
    geo_z->positions[0] = {0., 0., -1.};
    PointData<Point> ax_z_q(*ax_z);
    ax_z_q[0] = {0., 0., 0.25};
    auto pc_z = polyscope::registerPointCloud("World z-axis", geo_z->positions);
    pc_z->setPointColor({0., 0., 0.});
    auto vis_z = pc_z->addVectorQuantity("Axes", ax_z_q, VectorType::AMBIENT);
    vis_z->setVectorColor({0., 0., 1.});
    vis_z->setVectorRadius(0.01);
    vis_z->setEnabled(true);
}