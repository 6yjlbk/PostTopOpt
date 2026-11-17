#pragma once
#include "Mesh.h"
#include <functional>
namespace geometry {
struct Grid {
    Vec3 origin = Vec3(-3,-3,-3); // Начальная точка сетки
    Vec3 spacing = Vec3::Constant(0.1); // Шаг по каждой оси
    std::array<int,3> cells{80,80,80}; // Количество ячеек по каждой оси
};
using Field = std::function<double(const Vec3&)>; // Возвращает скалярное значение поля (заданное аналитически/через полигональную модель/как массив плотностей)
using Gradient = std::function<Vec3(const Vec3&)>; // Если она отсуствует, то градиент вычисляется численно 
void validateGrid(const Grid&);
std::size_t cellCount(const Grid&);
std::size_t cellIndex(const Grid&, int x, int y, int z); // z fastest
}
