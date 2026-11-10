#pragma once
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <vector>
namespace geometry {
using Vec3 = Eigen::Vector3d; // Вектор из 3х чисел типа вектор
struct Mesh {
    std::vector<Vec3> vertices; // Хранит координаты вершин поверхности
    std::vector<std::array<std::size_t, 4>> quads; // quads[] хранит четверки индексов этих вершин
    std::vector<std::array<std::size_t, 3>> triangles; // хранит список треугольников, каждый из которых задан тремя индексами вершин:
};
}
