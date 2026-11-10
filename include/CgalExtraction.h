#pragma once
#include "DualContouring.h"
namespace geometry {
struct CgalContourOptions {
    double isovalue = 0;
    bool constrainToCell = true;
    double gradientStep = 1e-5;
};
bool cgalAvailable() noexcept;
// CGAL uses its own vertex solver; no mapping of our QEF regularization.
Mesh cgalDualContour(const Grid&, const Field&, const Gradient& = {},
                     const CgalContourOptions& = {});
}
