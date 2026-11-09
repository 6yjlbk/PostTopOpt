#pragma once
#include "MarchingCubes.h"
#include <stdexcept>
#include <string_view>
namespace geometry {
enum class ExtractionMethod { DualContouring, EmilkDualContouring, MarchingCubes };
inline ExtractionMethod parseMethod(std::string_view method) {
    if (method == "dc") return ExtractionMethod::DualContouring;
    if (method == "emilk") return ExtractionMethod::EmilkDualContouring;
    if (method == "mc") return ExtractionMethod::MarchingCubes;
    throw std::invalid_argument("Unknown method; choose dc, emilk or mc");
}
inline Mesh extractSurface(ExtractionMethod method, const Grid& grid, const Field& field,
                           const Gradient& gradient = {}, const ContourOptions& options = {}) {
    switch (method) {
    case ExtractionMethod::DualContouring: return dualContour(grid, field, gradient, options);
    case ExtractionMethod::EmilkDualContouring: return dualContourEmilk(grid, field, gradient, options);
    case ExtractionMethod::MarchingCubes: return marchingCubes(grid, field, options.isovalue);
    }
    throw std::invalid_argument("Invalid extraction method");
}
}
