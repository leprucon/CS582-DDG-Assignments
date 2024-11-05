#include "include/ex2.h"
#include "include/utils.h"

bool near_zero(double x)
{
    return std::abs(x) < 1e-9;
}

double Exercise2::iso_value(const Point &_pt) const
{

    double iso = 0.;

    // TODO: ADD YOUR CODE HERE
    // 0 -> (x^2 + y^2)^0.5 - 1 = 0
    // 1 -> y^2 - sin(x^2) = 0
    // 2 -> sin(2x + 2y) - cos(4xy) + 1 = 0
    // 3 -> (3x^2 - y^2)^2 * y^2 - (x^2 + y^2)^4 = 0
    // 4 -> ? = 0

    if (function_id == 0)
    {
        iso = std::sqrt(std::pow(_pt.x, 2) + std::pow(_pt.y, 2)) - 1;
    }
    if (function_id == 1)
    {
        iso = std::pow(_pt.y, 2) - std::sin(std::pow(_pt.x, 2));
    }
    if (function_id == 2)
    {
        iso = std::sin(2 * _pt.x + 2 * _pt.y) - std::cos(4 * _pt.x * _pt.y) + 1;
    }
    if (function_id == 3)
    {
        iso = std::pow(3 * std::pow(_pt.x, 2) - std::pow(_pt.y, 2), 2) * std::pow(_pt.y, 2) - std::pow(std::pow(_pt.x, 2) + std::pow(_pt.y, 2), 4);
    }
    if (function_id == 4)
    {
        iso = std::pow(_pt.y - std::cbrt(std::pow(_pt.x, 2)), 2) + std::pow(_pt.x, 2) - 1;
    }
    return iso;
}

void Exercise2::compute_segment_points()
{
    segment_points.clear();

    // TODO: ADD YOUR CODE HERE
    // There is a check for segment_points to be of even length, so I assume we store redundant segment positions

    for (auto face : gc_mesh->faces())
    {
        geometry->requireVertexPositions();
        std::vector<Point> adjv;
        for (auto v : face.adjacentVertices())
        {
            adjv.push_back(geometry->vertexPositions[v]);
        }
        for (size_t i = 0; i < adjv.size(); i += 1)
        {
            auto value1 = iso_value(adjv[i]);
            auto value2 = iso_value(adjv[(i + 1) % adjv.size()]);
            if (value1 == 0 || near_zero(value1))
            {
                auto x1 = adjv[i].x;
                auto y1 = adjv[i].y;
                auto z1 = adjv[i].z;
                segment_points.push_back(Point{x1, y1, z1});
                segment_points.push_back(Point{x1, y1, z1});
            }
            if (value1 * value2 < 0)
            { // technically inefficient way to check but cleaner
                auto x1 = adjv[i].x;
                auto x2 = adjv[(i + 1) % adjv.size()].x;
                auto y1 = adjv[i].y;
                auto y2 = adjv[(i + 1) % adjv.size()].y;
                auto z1 = adjv[i].z;
                auto z2 = adjv[(i + 1) % adjv.size()].z;
                auto t = value1 / (value1 - value2);

                segment_points.push_back(Point{x1 * (1 - t) + x2 * t, y1 * (1 - t) + y2 * t, z1 * (1 - t) + z2 * t});
            }
        }
    }
}

void Exercise2::show_isovalue_and_level_set()
{
    show_isovalue();
    create_level_set0_segments();
}

void Exercise2::show_isovalue()
{

    std::vector<double> iso_values;
    for (auto v : gc_mesh->vertices())
    {
        Point point = geometry->vertexPositions[v];
        auto iv = iso_value(point);
        iso_values.push_back(iv);
    }

    auto max_iv = *std::max_element(iso_values.begin(), iso_values.end());
    auto min_iv = *std::min_element(iso_values.begin(), iso_values.end());

    const double range = max_iv - min_iv;
    std::vector<std::array<double, 3>> psColor;

    for (auto vh : gc_mesh->vertices())
    {
        auto t = (iso_values[vh.getIndex()] - min_iv) / range;
        auto color = util::map_val2color(t, 0, 1.0);
        psColor.push_back({color.x, color.y, color.z});
    }

    auto cq = ps_mesh->addVertexColorQuantity("Vertex Colors", psColor);
    cq->setEnabled(true);
}

void Exercise2::create_level_set0_segments()
{
    compute_segment_points();
    edges.clear();
    vertices.clear();

    if (segment_points.empty())
    {
        std::cerr << "segment_points is empty" << std::endl;
        return;
    }

    if (segment_points.size() % 2 != 0)
    {
        std::cerr << "Size of segment_points is not even" << std::endl;
        return;
    }

    for (size_t i = 0; i < segment_points.size(); i += 2)
    {
        auto p0 = segment_points[i];
        auto p1 = segment_points[i + 1];
        vertices.push_back({p0[0], p0[1], p0[2]});
        vertices.push_back({p1[0], p1[1], p1[2]});
        edges.push_back({i, i + 1});
    }

    polyscope::registerCurveNetwork("Level Set 0", vertices, edges);
}