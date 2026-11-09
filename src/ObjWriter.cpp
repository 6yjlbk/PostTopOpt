#include "ObjWriter.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace geometry {
void writeObj(const Mesh& mesh, const std::filesystem::path& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open OBJ output");
    out << std::setprecision(17);
    for (const auto& v : mesh.vertices) {
        if (!v.allFinite()) throw std::invalid_argument("Non-finite mesh vertex");
        out << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
    }
    for (const auto& q : mesh.quads) {
        out << "f";
        for (auto i : q) {
            if (i >= mesh.vertices.size()) throw std::invalid_argument("Invalid face index");
            out << ' ' << i + 1;
        }
        out << '\n';
    }
    for (const auto& triangle : mesh.triangles) {
        out << "f";
        for (auto i : triangle) {
            if (i >= mesh.vertices.size()) throw std::invalid_argument("Invalid face index");
            out << ' ' << i + 1;
        }
        out << '\n';
    }
    if (!out) throw std::runtime_error("Cannot write OBJ output");
}
}
