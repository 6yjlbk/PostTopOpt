#pragma once
#include "CellScalarVolume.h"
namespace geometry {
struct ParaViewContourOptions {
    double isovalue = 0.5; // Native cell values, no implicit 255 scaling.
    bool enableCapping = false;
    bool mergePoints = true;
};
bool paraViewAvailable() noexcept;
Mesh paraViewDualContour(const CellScalarVolume&, const ParaViewContourOptions& = {});
}
