#pragma once
#include "PlyImporter.h"
#include "DualContouring.h"
namespace geometry {
struct PlyDistanceField {
    Field field;
    Grid grid;
    PlyImportReport source;
};
// resolution: cells across the longest unpadded AABB side; padding: cells per side.
// Callbacks retain VTK objects. Use sequentially, not concurrently.
PlyDistanceField makePlyDistanceField(const std::filesystem::path&, int resolution = 64,
                                    int padding = 3);
}
