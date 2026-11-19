#pragma once
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace postop::vtkop_contract { // пространствно имен
struct UniformGridObject; // регулярная сетка
struct AmrObject; // контейнер сетки с VF, с которым работает алгоритм DC
struct ContourObject; // Алгоритм DC 
struct PolyDataObject; // Хранение результатов DC
using UniformGrid = std::shared_ptr<UniformGridObject>; 
using Amr = std::shared_ptr<AmrObject>;
using Contour = std::shared_ptr<ContourObject>;
using PolyData = std::shared_ptr<PolyDataObject>;

// Создание регулярной сетки
UniformGrid createUniformGrid(
    const std::array<double, 3>& origin, // начало сетки
    const std::array<double, 3>& spacing, // шаг сетки
    const std::array<int, 3>& cells, // количество ячеек по каждой оси
    std::size_t valueCount); // количество значений VF

// Заполнение ячеек значениями VF
void setCellDensities(const UniformGrid& grid,
                    const double* values, std::size_t count);

// Создания AMR-контейнера со значениями объемных долей
Amr createSingleBlockAmr(const UniformGrid& grid);

// Создание объекта фильтра DC
Contour createAmrDualContour();

// Настройка DC
void setDensityInput(const Contour& contour, const Amr& input); // Передача ячеек с виртуалной плотностью
void setIsovalue(const Contour& contour, double value); // Задание уровня Iso-поверхности
void setCapping(const Contour& contour, bool enabled); // Вкл/выкл закрытие поверхности на границе расчетной области
void setMergePoints(const Contour& contour, bool enabled); // Вкл/выкл объединение общих точек выходной сетки

// Запуск алгоритма
void update(const Contour& contour);

// Контейнер с результатами DC (список полигональных блоков результата работы алгоритма DC)
std::vector<PolyData> polyDataBlocks(const Contour& contour);

// Триангулирование результатов работы алгоритма (полигональная модель -> триангулирование представление)
PolyData triangulate(const PolyData& input);

// Чтение координат и граней
std::size_t pointCount(const PolyData& data); // Возвращает количество вершин новой поверхности
std::size_t triangleCount(const PolyData& data); // Возвращает количество граней вершин новой поверхности   
std::array<double, 3> point(const PolyData& data, std::size_t index); // Возвращает координаты заданой вершины
std::array<std::size_t, 3> triangle(const PolyData& data, std::size_t index); // Возвращает индексы трех узлов - вершин треугольника 
} // namespace postop::vtkop_contract
