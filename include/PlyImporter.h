#pragma once
#include "Mesh.h"
#include <filesystem>
namespace geometry {
struct PlyImportReport {
    std::size_t sourcePoints = 0;
    std::size_t sourcePolygons = 0;
    std::size_t boundaryEdges = 0;
    std::size_t nonManifoldEdges = 0;
    Vec3 minimum = Vec3::Zero();
    Vec3 maximum = Vec3::Zero();
    // Edge incidence only; not a proof of vertex manifoldness or absence of intersections.
    bool closedByEdgeIncidence() const { return boundaryEdges == 0 && nonManifoldEdges == 0; }
};
struct ImportedPly {
    Mesh mesh;
    PlyImportReport report;
};
// Requires the optional PostTopOptVtk library. Preserves coordinates and units.
// Does not merge points, fill holes, orient faces, normalize or build an SDF.
ImportedPly readPly(const std::filesystem::path& path);
}
