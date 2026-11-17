#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace postop {
struct RegularGrid {
    std::array<double, 3> origin{};
    std::array<double, 3> spacing{};
    std::array<int, 3> cells{};
};
struct ExtractionOptions {
    double threshold = 0.5;
    bool enableCapping = false;
    bool mergePoints = true;
};
struct Surface {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;
};
// Synchronous; input is borrowed only during this call, output owns its data.
// VF index: (x * Ny + y) * Nz + z. Values and threshold must be in [0,1].
// At least 3 cells per axis. Empty output is valid.
// Invalid input: std::invalid_argument; backend failure: std::runtime_error.
// Allocation failures propagate. No filesystem access.
Surface extractVolumeSurface(
    const RegularGrid& grid,
    const std::vector<double>& volumeFractions,
    const ExtractionOptions& options = {});
} // namespace postop
