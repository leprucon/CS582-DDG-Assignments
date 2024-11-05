#ifndef CURVESMOOTHING_HH
#define CURVESMOOTHING_HH

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include <chrono>

using namespace polyscope;
using geometrycentral::pointcloud::PointData;
using geometrycentral::pointcloud::PointPositionGeometry;

class Exercise3
{
public:
  using Point = glm::vec3;
  using PCPointCloud = geometrycentral::pointcloud::PointCloud;

  Exercise3() : Curve(nullptr) {}
  ~Exercise3() {}
  size_t num_vertices = 30;
  float epsilon = 0.01;
  int smooth_type = 0;
  int num_iter = 10;

  double maxX;
  double minX;
  double maxY;
  double minY;
  double maxZ;
  double minZ;
  bool first;
  void apply_smooth()
  {
    // std::cout << "Smoothing type: " << smooth_type << std::endl;
    if (smooth_type == 0)
    {
      laplacian_smoothing();
    }
    else if (smooth_type == 1)
    {
      osculating_circle();
    }
  };
  void generate_curve_simple();
  void generate_curve_figure_eight();
  void generate_curve_limacon();
  void generate_curve_3d();
  void generate_world_axes();

  void laplacian_smoothing();
  void osculating_circle();
  Point CircumscribedCircle(const Point &a, const Point &b, const Point &c);

  // taken from https://www.mcs.anl.gov/~fathom/meshkit-docs/html/circumcenter_8cpp_source.html
  void tricircumcenter3d(double a[3], double b[3], double c[3], double circumcenter[3])
  {
    double xba, yba, zba, xca, yca, zca;
    double balength, calength;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;

    /* Use coordinates relative to point `a' of the triangle. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba;
    calength = xca * xca + yca * yca + zca * zca;

    /* Cross product of these edges. */
#ifdef EXACT
    /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
    /*   to ensure a correctly signed (and reasonably accurate) result, */
    /*   avoiding any possibility of division by zero.                  */

    A[0] = b[1];
    A[1] = b[2];
    B[0] = c[1];
    B[1] = c[2];
    C[0] = a[1];
    C[1] = a[2];
    xcrossbc = orient2d(A, B, C);

    A[0] = c[0];
    A[1] = c[2];
    B[0] = b[0];
    B[1] = b[2];
    C[0] = a[0];
    C[1] = a[2];
    ycrossbc = orient2d(A, B, C);

    A[0] = b[0];
    A[1] = b[1];
    B[0] = c[0];
    B[1] = c[1];
    C[0] = a[0];
    C[1] = a[1];
    zcrossbc = orient2d(A, B, C);

    /*
    xcrossbc = orient2d(b[1], b[2], c[1], c[2], a[1], a[2]);
    ycrossbc = orient2d(b[2], b[0], c[2], c[0], a[2], a[0]);
    zcrossbc = orient2d(b[0], b[1], c[0], c[1], a[0], a[1]);
    */
#else
    printf(" Warning: IEEE floating points used: Define -DEXACT in makefile \n");
    /* Take your chances with floating-point roundoff. */
    xcrossbc = yba * zca - yca * zba;
    ycrossbc = zba * xca - zca * xba;
    zcrossbc = xba * yca - xca * yba;
#endif

    /* Calculate the denominator of the formulae. */
    denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
                         zcrossbc * zcrossbc);

    /* Calculate offset (from `a') of circumcenter. */
    xcirca = ((balength * yca - calength * yba) * zcrossbc -
              (balength * zca - calength * zba) * ycrossbc) *
             denominator;
    ycirca = ((balength * zca - calength * zba) * xcrossbc -
              (balength * xca - calength * xba) * zcrossbc) *
             denominator;
    zcirca = ((balength * xca - calength * xba) * ycrossbc -
              (balength * yca - calength * yba) * xcrossbc) *
             denominator;
    circumcenter[0] = xcirca;
    circumcenter[1] = ycirca;
    circumcenter[2] = zcirca;
  }

private:
  CurveNetwork *Curve;
  std::vector<Point> points;
  std::vector<std::array<size_t, 2>> edges;
};

#endif
