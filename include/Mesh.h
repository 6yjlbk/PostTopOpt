#pragma once
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <vector>
namespace geometry {
using Vec3 = Eigen::Vector3d;
struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<std::size_t, 4>> quads;
};
}
