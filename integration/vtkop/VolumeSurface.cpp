#include "VolumeSurface.h"
#include "VtkopContract.h"
#include <stdexcept>

namespace {
Surface extractAmrSurface(const Amr& amr,
                          const ExtractionOptions& options) {
    // Создание фильтра DC (создание объекта vtkAMRDualContour)
    auto contour = createAmrDualContour();
    setDensityInput(contour, amr); // Подключение поисковой сетки
    setIsovalue(contour, options.threshold); // Задание уровня изоповерхности
    setCapping(contour, options.enableCapping); // Настраиваем закрытие поверхности на границе области
    setMergePoints(contour, options.mergePoints); // Настраиваем объединение точек, предусмотренное фильтром
    update(contour, options.extractionSeconds); // Выполнение алгоритма DC

    Surface result; // Объект результатов
    for (const auto& poly : polyDataBlocks(contour)) { // Обход по полигональным блокам и запись его данных в result
        auto data = triangulate(poly);
        const auto points_prev = result.vertices.size(); //предыдущее количество вершин, нужно для корректного присваения вершин (актуально, если используется несколько блоков в организации данных VTK)
        const auto points = pointCount(data);
        const auto faces = triangleCount(data);
        if (points > result.vertices.max_size() - points_prev ||
            faces > result.triangles.max_size() - result.triangles.size())
            throw std::runtime_error("Surface exceeds container limits");
        result.vertices.reserve(points_prev + points);
        result.triangles.reserve(result.triangles.size() + faces);
        for (std::size_t i = 0; i < points; ++i) {
            const auto p = point(data, i);
            result.vertices.push_back(p);
        }
        for (std::size_t i = 0; i < faces; ++i) {
            auto face = triangle(data, i);
            for (auto& id : face) {
                id += points_prev;
            }
            result.triangles.push_back(face);
        }
    }
    return result;
}
} // namespace

Surface extractVolumeSurface(const std::vector<std::array<double, 3>>& points,
                             const std::vector<std::size_t>& offsets,
                             const std::vector<std::size_t>& nodes,
                             const std::vector<double>& volumeFractions,
                             const std::array<int, 3>& searchCells,
                             const ExtractionOptions& options) {
    auto amr = transferCellDensities( // Создание поисковой сетки и заполнение ее ячеек значениями VF
        points, offsets, nodes, volumeFractions, searchCells, options.transferMethod);
    return extractAmrSurface(amr, options);
}
