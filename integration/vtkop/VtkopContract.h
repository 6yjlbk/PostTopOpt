#pragma once
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

enum class TransferMethod { CellCenterSampling, NodalInterpolation, VolumeWeighted };

struct AmrObject; // контейнер сетки с VF, с которым работает алгоритм DC
struct ContourObject; // Фильтр DC
struct PolyDataObject; // Хранение результатов DC
using Amr = std::shared_ptr<AmrObject>;
using Contour = std::shared_ptr<ContourObject>;
using PolyData = std::shared_ptr<PolyDataObject>;


// Создание поисковой сетки и заполнение ее ячеек значениями VF
Amr transferCellDensities(
    const std::vector<std::array<double, 3>>& points,
    const std::vector<std::size_t>& offsets,
    const std::vector<std::size_t>& nodes,
    const std::vector<double>& values,
    const std::array<int, 3>& searchCells,
    TransferMethod method);

// Создание объекта фильтра DC
Contour createAmrDualContour();

// Настройка фильтра DC
void setDensityInput(const Contour& contour, const Amr& input); // Передача ячеек с виртуалной плотностью
void setIsovalue(const Contour& contour, double value); // Задание уровня Iso-поверхности
void setCapping(const Contour& contour, bool enabled); // Вкл/выкл закрытие поверхности на границе расчетной области
void setMergePoints(const Contour& contour, bool enabled); // Вкл/выкл объединение общих точек выходной сетки

// Запуск алгоритма
void update(const Contour& contour, double* extractionSeconds = nullptr);

// Контейнер с результатами DC (список полигональных блоков результата работы алгоритма DC)
std::vector<PolyData> polyDataBlocks(const Contour& contour);

// Триангулирование результатов работы алгоритма (полигональная модель -> триангулирование представление)
PolyData triangulate(const PolyData& input);

// Чтение координат и граней
std::size_t pointCount(const PolyData& data); // Возвращает количество вершин новой поверхности
std::size_t triangleCount(const PolyData& data); // Возвращает количество граней вершин новой поверхности
std::array<double, 3> point(const PolyData& data, std::size_t index); // Возвращает координаты заданой вершины
std::array<std::size_t, 3> triangle(const PolyData& data, std::size_t index); // Возвращает индексы трех узлов - вершин треугольника

#include <string>
#include <vtkSmartPointer.h>

class vtkUnstructuredGrid;
class vtkNonOverlappingAMR;
// Объявление объекта поисковой сетки
struct TransferGrid {
    std::array<double,3> origin{};
    std::array<double,3> spacing{};
    std::array<int,3> cells{};
};
// Получение поисковой сетки из исходной
TransferGrid makeTransferGrid(vtkUnstructuredGrid* source,
    const std::array<int,3>& cells);
// Перенос VF на поисковую сетку
vtkSmartPointer<vtkNonOverlappingAMR> transferVF(vtkUnstructuredGrid* source,
    const std::string& vfArray, const TransferGrid& grid, TransferMethod method); // source - исходная КЭ сетка в формате VTK, vfArray - массив исходных значений vf, grid - поисковая сетка, method - выбранный метод переноса
