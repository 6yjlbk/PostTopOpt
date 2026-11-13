#pragma once
#include "Field.h"
namespace geometry {
struct ContourOptions {
    double isovalue = 0;
    double regularization = 0.001;
    bool constrainToCell = true;
    double gradientStep = 1e-5;
    bool bisection = true;
    bool regularizationEnabled = true;
};
Mesh dualContour(const Grid&, const Field&, const Gradient& = {}, const ContourOptions& = {});
Mesh dualContourEmilk(const Grid&, const Field&, const Gradient& = {}, const ContourOptions& = {});
}
