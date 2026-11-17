#include "VolumeSurface.h"
#include "VtkopContract.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace postop {
namespace {
std::size_t validate(const RegularGrid& grid,
                     const std::vector<double>& values,
                     const ExtractionOptions& options) {
    std::size_t count = 1;
    for (std::size_t d = 0; d < 3; ++d) {
        const int n = grid.cells[d];
        if (n < 3 || n == (std::numeric_limits<int>::max)())
            throw std::invalid_argument("Expected 3..INT_MAX-1 cells per axis");
        const auto size = static_cast<std::size_t>(n);
        if (count > (std::numeric_limits<std::size_t>::max)() / size)
            throw std::invalid_argument("Grid cell count overflow");
        count *= size;
        if (!std::isfinite(grid.origin[d]) ||
            !std::isfinite(grid.spacing[d]) || grid.spacing[d] <= 0)
            throw std::invalid_argument("Invalid grid origin or spacing");
        const double upper = grid.origin[d] + grid.spacing[d] * n;
        if (!std::isfinite(upper) || upper <= grid.origin[d])
            throw std::invalid_argument("Invalid grid bounds");
    }
    if (values.size() != count)
        throw std::invalid_argument("VF count does not match the grid");
    for (double value : values)
        if (!std::isfinite(value) || value < 0 || value > 1)
            throw std::invalid_argument("VF must be finite and in [0,1]");
    if (!std::isfinite(options.threshold) ||
        options.threshold < 0 || options.threshold > 1)
        throw std::invalid_argument("Threshold must be finite and in [0,1]");
    return count;
}

template<class T>
void requireHandle(const T& handle) {
    if (!handle) throw std::runtime_error("vtkop returned a null object");
}
} // namespace

Surface extractVolumeSurface(const RegularGrid& grid,
                             const std::vector<double>& volumeFractions,
                             const ExtractionOptions& options) {
    namespace api = vtkop_contract;
    const auto count = validate(grid, volumeFractions, options);
    auto block = api::createUniformGrid(grid.origin, grid.spacing,
        {grid.cells[0] + 1, grid.cells[1] + 1, grid.cells[2] + 1});
    requireHandle(block);
    {
        // Retain the current z-fastest -> VTK x-fastest mapping.
        const auto nx = static_cast<std::size_t>(grid.cells[0]);
        const auto ny = static_cast<std::size_t>(grid.cells[1]);
        const auto nz = static_cast<std::size_t>(grid.cells[2]);
        std::vector<double> vtkValues(count);
        for (std::size_t z = 0; z < nz; ++z)
            for (std::size_t y = 0; y < ny; ++y)
                for (std::size_t x = 0; x < nx; ++x)
                    vtkValues[(z * ny + y) * nx + x] =
                        volumeFractions[(x * ny + y) * nz + z];
        api::setCellScalars(block, "density", vtkValues.data(), count);
    }
    auto amr = api::createSingleBlockAmr(block);
    requireHandle(amr);
    auto contour = api::createAmrDualContour();
    requireHandle(contour);
    api::setInput(contour, amr);
    api::setCellArray(contour, "density");
    api::setIsovalue(contour, options.threshold);
    api::setCapping(contour, options.enableCapping);
    api::setMergePoints(contour, options.mergePoints);
    api::update(contour);

    Surface result;
    for (const auto& poly : api::polyDataBlocks(contour)) {
        requireHandle(poly);
        auto data = api::triangulate(poly);
        requireHandle(data);
        const auto offset = result.vertices.size();
        const auto points = api::pointCount(data);
        const auto faces = api::triangleCount(data);
        if (points > result.vertices.max_size() - offset ||
            faces > result.triangles.max_size() - result.triangles.size())
            throw std::runtime_error("Surface exceeds container limits");
        result.vertices.reserve(offset + points);
        result.triangles.reserve(result.triangles.size() + faces);
        for (std::size_t i = 0; i < points; ++i) {
            const auto p = api::point(data, i);
            for (double coordinate : p)
                if (!std::isfinite(coordinate))
                    throw std::runtime_error("Non-finite surface coordinate");
            result.vertices.push_back(p);
        }
        for (std::size_t i = 0; i < faces; ++i) {
            auto face = api::triangle(data, i);
            for (auto& id : face) {
                if (id >= points)
                    throw std::runtime_error("Invalid surface vertex index");
                id += offset;
            }
            result.triangles.push_back(face);
        }
    }
    return result;
}
} // namespace postop
