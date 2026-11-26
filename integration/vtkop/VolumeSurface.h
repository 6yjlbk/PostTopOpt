#pragma once
#include <array>
#include <cstddef>
#include <vector>
#include "VtkopContract.h"


// Параметры работы алгоритма
struct ExtractionOptions {
    double iso = 0.5; // уровень изоповерхности
    bool enableCapping = true; // включить закрытие поверхности на границе расчётной области;
    bool mergePoints = true; // включить предусмотренное фильтром объединение точек
    TransferMethod transferMethod = TransferMethod::CellCenterSampling; // выбор метода
    double* extractionSeconds = nullptr; // Optional output: vtkAMRDualContour::Update only, seconds.
};

// Под замену пойдет, тип сделать изоповерхность
struct Surface {
    std::vector<std::array<double, 3>> vertices; //координаты вершин полученной поверхности
    std::vector<std::array<std::size_t, 3>> triangles; //тройки индексов в vertices, например triangles[0] = {2, 5, 7}
};
// Функция запуска алгоритма DC
Surface extractVolumeSurface(
    const std::vector<std::array<double, 3>>& points, // Координаты узлов исходной КЭ-сетки
    const std::vector<std::size_t>& offsets, // offsets отмечает, где начинается и заканчивается каждый элемент
    const std::vector<std::size_t>& nodes, // Индексы узлов всех элементов подряд
    const std::vector<double>& volumeFractions, // Значения VF исходных элементов в том же порядке
    const std::array<int, 3>& searchCells, // Размерность поисковой сетки (число ячеек по х, число ячеек по у, число ячеек по z). Границы AABB определяются по крайним координатам тела
    const ExtractionOptions& options = {}); // сейчас значения по умолчанию
