#pragma once
#include "Mesh.h"
#include <functional>
namespace geometry {
struct Grid {
    Vec3 origin = Vec3(-3, -3, -3); // Начало СК, из которой строится AABB
    Vec3 spacing = Vec3::Constant(0.1); // Шаг сетки
    std::array<int, 3> cells{80, 80, 80}; //Кол-во элементов по каждой из осей
};
struct ContourOptions {
    double isovalue = 0; // Уровень поля, поверхность которого извлекается f(x) = isovalue; далее f(x) вычитается из значений поля в каждой ячейки
    double regularization = 1e-6; // Штраф при размещении вершины вдали от ее среднеого положения
    bool constrainToCell = true; // Требование оставить вершину внутри своей ячейки
    double gradientStep = 1e-5; // 	Шаг численного вычисления градиента, если он не передан
};
using Field = std::function<double(const Vec3&)>; // Функция скалярного поля
using Gradient = std::function<Vec3(const Vec3&)>; // Функция градиента
// Negative values relative to isovalue denote the interior. Boundaries are not capped.
Mesh dualContour(const Grid&, const Field&, const Gradient& = {}, const ContourOptions& = {});
}
