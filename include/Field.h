#pragma once
#include <Eigen/Core>
#include <array>
#include <vector>
#include <functional>
#include <limits>
#include <stdexcept>
namespace geometry {
using Vec3 = Eigen::Vector3d; // Вектор из 3х чисел типа вектор
struct Mesh {
    std::vector<Vec3> vertices; // Хранит координаты вершин поверхности
    std::vector<std::array<std::size_t, 4>> quads; // quads[] хранит четверки индексов этих вершин
    std::vector<std::array<std::size_t, 3>> triangles; // хранит список треугольников, каждый из которых задан тремя индексами вершин
};
}

namespace geometry {
struct Grid {
    Vec3 origin = Vec3(-3,-3,-3); // Начальная точка сетки
    Vec3 spacing = Vec3::Constant(0.1); // Шаг по каждой оси
    std::array<int,3> cells{80,80,80}; // Количество ячеек по каждой оси
};
using Field = std::function<double(const Vec3&)>; // Возвращает скалярное значение поля (заданное аналитически/через полигональную модель/как массив плотностей)
using Gradient = std::function<Vec3(const Vec3&)>; // Если она отсуствует, то градиент вычисляется численно 
inline void validateGrid(const Grid&);
inline std::size_t cellCount(const Grid&);
inline std::size_t cellIndex(const Grid&, int x, int y, int z); // z fastest
}

namespace geometry {
inline std::size_t cellCount(const Grid& g) {
    std::size_t count=1;
    for(int n:g.cells) {
        if(n<=0 || n==(std::numeric_limits<int>::max)() || count>(std::numeric_limits<std::size_t>::max)()/static_cast<std::size_t>(n))
            throw std::invalid_argument("Invalid grid dimensions");
        count*=static_cast<std::size_t>(n);
    }
    return count;
}
inline void validateGrid(const Grid& g) {
    cellCount(g);
    if(!g.origin.allFinite() || !g.spacing.allFinite() || (g.spacing.array()<=0).any())
        throw std::invalid_argument("Invalid grid geometry");
    const Vec3 upper=g.origin+g.spacing.cwiseProduct(Vec3(g.cells[0],g.cells[1],g.cells[2]));
    if(!upper.allFinite() || (upper.array()<=g.origin.array()).any()) throw std::invalid_argument("Invalid grid bounds");
}
inline std::size_t cellIndex(const Grid& g,int x,int y,int z) {
    if(x<0 || y<0 || z<0 || x>=g.cells[0] || y>=g.cells[1] || z>=g.cells[2]) throw std::out_of_range("Cell index");
    return (static_cast<std::size_t>(x)*g.cells[1]+y)*g.cells[2]+z;
}
}



