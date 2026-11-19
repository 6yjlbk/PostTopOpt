#include "VolumeSurface.h"
#include "VtkopContract.h"
#include <stdexcept>

namespace postop {
Surface extractVolumeSurface(const RegularGrid& grid,
                             const std::vector<double>& volumeFractions,
                             const ExtractionOptions& options) {
    namespace api = vtkop_contract;
    // Создание сетки
    const auto count = volumeFractions.size();
    auto block = api::createUniformGrid(grid.origin, grid.spacing, grid.cells, count);
    // блок с временными переменными, которые уничтожаются в конце блока
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
        api::setCellDensities(block, vtkValues.data(), count);
    }
    auto amr = api::createSingleBlockAmr(block);
    auto contour = api::createAmrDualContour();
    api::setDensityInput(contour, amr);
    api::setIsovalue(contour, options.threshold);
    api::setCapping(contour, options.enableCapping);
    api::setMergePoints(contour, options.mergePoints);
    api::update(contour);

    Surface result;
    for (const auto& poly : api::polyDataBlocks(contour)) {
        auto data = api::triangulate(poly);
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
            result.vertices.push_back(p);
        }
        for (std::size_t i = 0; i < faces; ++i) {
            auto face = api::triangle(data, i);
            for (auto& id : face) {
                id += offset;
            }
            result.triangles.push_back(face);
        }
    }
    return result;
}
} // namespace postop
