#pragma once
#include "DualContouring.h"
namespace geometry {
// Classic 256-case Marching Cubes. No QEF or gradient queries.
// Ambiguous cases use the fixed table, not MC33/asymptotic-decider topology.
Mesh marchingCubes(const Grid&, const Field&, double isovalue = 0);
}
