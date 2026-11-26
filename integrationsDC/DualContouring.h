#pragma once
#include "BuildIsoSurface.h"

namespace postopdc {
// Temporary FE field implementation, ready to move into the host element class.
// q is a NATURAL coordinate, not a physical position:
// tetra: r,s,t >= 0, r+s+t <= 1; wedge: r,s >= 0, r+s <= 1, 0 <= t <= 1;
// hexa: [0,1]^3. gradient(q) returns derivatives in PHYSICAL coordinates.
class ElementField {
public:
    ElementField(const Element* element, const std::vector<Vec3>& points,
                 const std::vector<float>& nodalVF)
        : element(element), points(points), nodalVF(nodalVF) {}
    Vec3 position(const Vec3& q) const;
    float value(const Vec3& q) const;
    Vec3 gradient(const Vec3& q) const;
    Eigen::Matrix3f jacobian(const Vec3& q) const;
    Vec3 project(const Vec3& q) const;
    const Element* element;
private:
    void shape(const Vec3& q, std::array<float,8>& n,
               std::array<Vec3,8>& derivatives) const;
    const std::vector<Vec3>& points;
    const std::vector<float>& nodalVF;
};
Surface dualContour(const std::vector<Vec3>& points,
    const std::vector<const Element*>& elements,
    const std::vector<float>& nodalVF, const ExtractionOptions& options);
} // namespace postopdc
