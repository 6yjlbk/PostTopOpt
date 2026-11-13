#pragma once
#include "CellScalarVolume.h"
namespace geometry {
struct ParaViewProfile { double preparationMs=0, filterMs=0, conversionMs=0; };
struct ParaViewContourOptions {
    double isovalue = 0.5; // Native cell values, no implicit 255 scaling.
    bool enableCapping = false;
    bool mergePoints = true;
    ParaViewProfile* profile=nullptr;
};
bool paraViewAvailable() noexcept;
Mesh paraViewDualContour(const CellScalarVolume&, const ParaViewContourOptions& = {});
}
