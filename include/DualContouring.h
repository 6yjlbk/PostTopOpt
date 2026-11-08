#pragma once
#include "Mesh.h"
#include <functional>
namespace geometry {
struct Grid {
    Vec3 origin = Vec3(-2, -2, -2);
    Vec3 spacing = Vec3::Constant(0.1);
    std::array<int, 3> cells{40, 40, 40};
};
struct ContourOptions {
    double isovalue = 0;
    double regularization = 1e-6;
    bool constrainToCell = true;
    double gradientStep = 1e-5;
};
using Field = std::function<double(const Vec3&)>;
using Gradient = std::function<Vec3(const Vec3&)>;
// Negative values relative to isovalue denote the interior. Boundaries are not capped.
Mesh dualContour(const Grid&, const Field&, const Gradient& = {}, const ContourOptions& = {});
}
