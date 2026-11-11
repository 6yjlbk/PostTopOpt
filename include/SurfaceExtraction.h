#pragma once
#include "MarchingCubes.h"
#include "CgalExtraction.h"
#include <stdexcept>
#include <string_view>
namespace geometry {
enum class ExtractionMethod { DualContouring, EmilkDualContouring, MarchingCubes, CgalDualContouring, CgalOctreeDualContouring };
inline ExtractionMethod parseMethod(std::string_view method) {
    if (method == "dc") return ExtractionMethod::DualContouring;
    if (method == "emilk") return ExtractionMethod::EmilkDualContouring;
    if (method == "cgal-dc-octree") return ExtractionMethod::CgalOctreeDualContouring;
    if (method == "cgal-dc") return ExtractionMethod::CgalDualContouring;
    if (method == "mc") return ExtractionMethod::MarchingCubes;
    throw std::invalid_argument("Unknown method; choose dc, emilk, mc, cgal-dc or cgal-dc-octree");
}
inline Mesh extractSurface(ExtractionMethod method, const Grid& grid, const Field& field,
                           const Gradient& gradient = {}, const ContourOptions& options = {}) { 
    switch (method) { //РІС‹Р·РѕРІ РЅСѓР¶РЅРѕР№ С„СѓРЅРєС†РёРё
    case ExtractionMethod::DualContouring: return dualContour(grid, field, gradient, options);
    case ExtractionMethod::EmilkDualContouring: return dualContourEmilk(grid, field, gradient, options);
    case ExtractionMethod::CgalDualContouring:
        // CGAL has its own vertex solver; ContourOptions::regularization is not used.
        return cgalDualContour(grid, field, gradient,
            {options.isovalue, options.constrainToCell, options.gradientStep});
    case ExtractionMethod::CgalOctreeDualContouring: {
        CgalContourOptions cgalOptions;
        cgalOptions.isovalue=options.isovalue;
        cgalOptions.constrainToCell=options.constrainToCell;
        cgalOptions.gradientStep=options.gradientStep;
        cgalOptions.octree=true;
        return cgalDualContour(grid,field,gradient,cgalOptions);
    }
    case ExtractionMethod::MarchingCubes: return marchingCubes(grid, field, options.isovalue);
    }
    throw std::invalid_argument("Invalid extraction method");
}
}
