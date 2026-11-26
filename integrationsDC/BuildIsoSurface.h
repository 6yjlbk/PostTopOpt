#pragma once
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <vector>

namespace postopdc {
using Vec3 = Eigen::Vector3f;
enum class ElementType { Tetra4 = 4, Wedge6 = 6, Hexa8 = 8 };
// Linear conforming FE mesh; node order follows the conventional VTK ordering.
// Only the first nodeCount() entries are used. Geometry remains owned by the caller.
struct Element {
    ElementType type = ElementType::Tetra4;
    std::array<std::size_t,8> nodes{};
    float vf = 0;
    int nodeCount() const { return static_cast<int>(type); }
};
struct ExtractionOptions {
    float iso = 0.5f;
    float regularization = 0.001f;
    bool regularizationEnabled = true;
    bool constrainToCell = true;
    bool bisection = true;
    double* extractionSeconds = nullptr; // DC only, excludes nodal averaging.
};
struct Surface {
    std::vector<std::array<float,3>> vertices;
    std::vector<std::array<std::size_t,4>> quads;
    std::vector<std::array<std::size_t,3>> triangles;
    // General edge valence needs polygons; do not silently fan-triangulate them.
    std::vector<std::vector<std::size_t>> polygons;
};
// No padding or boundary caps. All incident elements must be supplied.
// Pointers are borrowed for this call; input elements are not modified/copied.
Surface buildIsoSurface(const std::vector<Vec3>& points,
    const std::vector<const Element*>& elements,
    const ExtractionOptions& options = {});
} // namespace postopdc
