#include <cctype>
#include <string>
#include <stdexcept>
#include "PlyImporter.h"
#include "ObjWriter.h"
#include <algorithm>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    try {
        std::vector<std::filesystem::path> files;
        std::filesystem::path output;
        for (int i=1;i<argc;++i) {
            std::string arg=argv[i];
            if (arg=="--help") {
                std::cout << "Usage: PostTopOptPlyImport [file.ply ...] [--output-dir directory]\n"
                             "Without files: inspect top-level examples/*.ply only. No SDF or DC run.\n";
                return 0;
            }
            if (arg=="--output-dir") {
                if (++i==argc) throw std::invalid_argument("Missing output directory");
                output=argv[i];
            } else if (arg.starts_with("--")) throw std::invalid_argument("Unknown option: "+arg);
            else files.emplace_back(arg);
        }
        if (files.empty()) {
            for (const auto& entry:std::filesystem::directory_iterator(POSTTOPOPT_EXAMPLES_DIR)) {
                auto extension=entry.path().extension().string();
                std::transform(extension.begin(),extension.end(),extension.begin(),
                    [](unsigned char c){return static_cast<char>(std::tolower(c));});
                if (entry.is_regular_file() && extension==".ply") files.push_back(entry.path());
            }
        }
        std::sort(files.begin(),files.end());
        if (files.empty()) throw std::runtime_error("No PLY files found");
        bool failed=false;
        for (const auto& file:files) {
            try {
                const auto imported=geometry::readPly(file);
                const auto& r=imported.report;
                std::cout << file.string() << "\nsource_points=" << r.sourcePoints
                    << ", source_polygons=" << r.sourcePolygons
                    << ", triangles=" << imported.mesh.triangles.size()
                    << ", boundary_edges=" << r.boundaryEdges
                    << ", non_manifold_edges=" << r.nonManifoldEdges
                    << ", closed_by_edge_incidence=" << std::boolalpha << r.closedByEdgeIncidence()
                    << "\nbounds_min=" << r.minimum.transpose() << ", bounds_max=" << r.maximum.transpose() << '\n';
                if (!output.empty()) {
                    std::filesystem::create_directories(output);
                    const auto destination=output/(file.stem().string()+"_import.obj");
                    if (std::filesystem::exists(destination)) throw std::runtime_error("Output exists: "+destination.string());
                    geometry::writeObj(imported.mesh,destination);
                }
            } catch (const std::exception& e) { failed=true; std::cerr<<file.string()<<": "<<e.what()<<'\n'; }
        }
        return failed ? 1 : 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
