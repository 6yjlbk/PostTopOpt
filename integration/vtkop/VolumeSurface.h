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

// Под замену пойдет, тип сделать изоповерхность
struct Surface {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<std::size_t, 3>> triangles;
};
// Функция, внутри которой происходит запуск алгоритма и вся работа над ним
Surface extractVolumeSurface(
    const RegularGrid& grid,
    const std::vector<double>& volumeFractions,
    const ExtractionOptions& options = {});
} // namespace postop
