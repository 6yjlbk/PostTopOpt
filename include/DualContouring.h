#pragma once
#include "Field.h"
namespace geometry {
struct DcProfile {
    double samplingMs=0, hermiteMs=0, qefMs=0, connectivityMs=0;
    std::size_t fieldQueries=0, activeCells=0, edgeVisits=0, uniqueEdges=0;
    std::size_t hermiteEvaluations=0, hermiteCacheHits=0;
    std::size_t qefInteriorReturns=0, qefBoundaryCells=0;
};
struct ContourOptions {
    double isovalue = 0;
    double regularization = 0.001;
    bool constrainToCell = true;
    double gradientStep = 1e-5;
    bool bisection = true;
    bool regularizationEnabled = true;
    bool optimizeQef = true; // Control path remains available for numerical comparisons.
    bool cacheHermite = true; // Deterministic Field/Gradient required; classic DC only.
    DcProfile* profile = nullptr; // Optional diagnostic run; reset on entry.
};
Mesh dualContour(const Grid&, const Field&, const Gradient& = {}, const ContourOptions& = {});
Mesh dualContourEmilk(const Grid&, const Field&, const Gradient& = {}, const ContourOptions& = {});
}
