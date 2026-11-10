#pragma once
#include "SurfaceExtraction.h"
#include "ObjWriter.h"
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>

namespace examples {
inline geometry::Grid grid(double halfExtent, int cells) {
    geometry::Grid result;
    result.origin = geometry::Vec3::Constant(-halfExtent); // Начало СК, из которой строится AABB
    result.spacing = geometry::Vec3::Constant(2 * halfExtent / cells); // Шаг сетки
    result.cells = {cells, cells, cells}; //Кол-во элементов по каждой из осей
    return result;
}

inline int run(int argc, char** argv, const char* defaultPath,
               const geometry::Grid& grid, const geometry::Field& field,
               const geometry::Gradient& gradient = {}) {
    try {
        std::string methodName="dc";
        std::filesystem::path path;
        for (int i=1; i<argc; ++i) {
            const std::string argument=argv[i];
            if (argument=="--help") {
                std::cout << "Usage: executable [output.obj] [--method dc|emilk|mc|cgal-dc]\n";
                return 0;
            }
            if (argument=="--method") {
                if (++i==argc) throw std::invalid_argument("--method requires dc, emilk, mc or cgal-dc");
                methodName=argv[i];
            } else {
                if (argument.starts_with("--") || !path.empty())
                    throw std::invalid_argument("Usage: executable [output.obj] [--method dc|emilk|mc|cgal-dc]");
                path=argument;
            }
        }
        const auto method=geometry::parseMethod(methodName);
        if (path.empty()) {
            path=std::filesystem::path(POSTTOPOPT_RESULTS_DIR)/defaultPath;
            // Preserve the old DC default; separate the new methods' outputs.
            if (methodName!="dc")
                path=path.parent_path()/(path.stem().string()+"_"+methodName+path.extension().string());
        }
        const auto begin=std::chrono::steady_clock::now();
        const auto mesh=geometry::extractSurface(method,grid,field,gradient);
        const auto elapsed=std::chrono::duration<double,std::milli>(
            std::chrono::steady_clock::now()-begin).count();
        if (mesh.vertices.empty() || (mesh.quads.empty() && mesh.triangles.empty()))
            throw std::runtime_error("Example produced an empty surface");
        if (!path.parent_path().empty())
            std::filesystem::create_directories(path.parent_path());
        geometry::writeObj(mesh,path);
        std::cout << methodName << ": " << mesh.vertices.size() << " vertices, "
                  << mesh.quads.size() << " quads, " << mesh.triangles.size()
                  << " triangles, extraction_ms=" << elapsed << " -> " << path.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
}
