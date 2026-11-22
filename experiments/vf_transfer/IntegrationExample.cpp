#include "VolumeSurface.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// This calling application owns its FE data. Only VolumeSurface.h is used.
int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3) throw std::runtime_error("Usage: IntegrationExample <new-output-directory> [shape]");
        const std::filesystem::path output = argv[1];
        if (std::filesystem::exists(output)) throw std::runtime_error("Output directory already exists");
        std::filesystem::create_directories(output);
        constexpr int sourceCells = 50;
        constexpr double step = 2.0 / sourceCells;
        std::vector<std::array<double, 3>> points;
        for (int z = 0; z <= sourceCells; ++z)
            for (int y = 0; y <= sourceCells; ++y)
                for (int x = 0; x <= sourceCells; ++x) {
                    std::array<double, 3> point{-1 + x * step, -1 + y * step, -1 + z * step};
                    if (x > 0 && x < sourceCells && y > 0 && y < sourceCells && z > 0 && z < sourceCells)
                        for (int axis = 0; axis < 3; ++axis)
                            point[axis] += .06 * step * std::sin(1.7*x + 2.3*y + 3.1*z + 1.2*axis);
                    points.push_back(point);
                }
        const auto pointId = [](int x, int y, int z) -> std::size_t {
            return (std::size_t(z) * (sourceCells + 1) + y) * (sourceCells + 1) + x;
        };
        std::vector<std::size_t> nodes;
        std::vector<std::size_t> offsets{0};
        std::vector<std::array<double, 3>> centers;
        for (int z = 0; z < sourceCells; ++z)
            for (int y = 0; y < sourceCells; ++y)
                for (int x = 0; x < sourceCells; ++x) {
                    std::array<int, 3> permutation{0, 1, 2};
                    do {
                        int position[3]{x, y, z};
                        std::array<std::size_t, 4> ids;
                        ids[0] = pointId(x, y, z);
                        for (int axis = 0; axis < 3; ++axis) {
                            ++position[permutation[axis]];
                            ids[axis + 1] = pointId(position[0], position[1], position[2]);
                        }
                        // Alternate permutations reverse orientation.
                        if ((permutation[0] > permutation[1]) + (permutation[0] > permutation[2]) +
                            (permutation[1] > permutation[2]) & 1)
                            std::swap(ids[1], ids[2]);
                        std::array<double, 3> center{};
                        for (auto id : ids)
                            for (int axis = 0; axis < 3; ++axis) center[axis] += points[id][axis] * .25;
                        nodes.insert(nodes.end(), ids.begin(), ids.end());
                        offsets.push_back(nodes.size());
                        centers.push_back(center);
                    } while (std::next_permutation(permutation.begin(), permutation.end()));
                }
        std::ofstream report(output / "runs.csv");
        report << "shape,method,vertices,triangles,extractDC_seconds,extractVolumeSurface_seconds\n";
        for (const std::string shape : {"sphere", "cylinder", "torus", "gyroid"}) {
            if (argc == 3 && shape != argv[2]) continue;
            std::vector<double> fractions;
            for (const auto& center : centers) {
                const auto [x, y, z] = center;
                double distance;
                if (shape == "sphere") distance = std::sqrt(x*x + y*y + z*z) - .65;
                else if (shape == "cylinder") distance = std::max(std::hypot(x, y) - .5, std::abs(z) - .65);
                else if (shape == "torus") distance = std::hypot(std::hypot(x, y) - .48, z) - .18;
                else {
                    constexpr double k = 6.283185307179586;
                    distance = std::max((std::abs(std::sin(k*x)*std::cos(k*y) + std::sin(k*y)*std::cos(k*z) +
                        std::sin(k*z)*std::cos(k*x)) - .7) / k,
                        std::max({std::abs(x), std::abs(y), std::abs(z)}) - .78);
                }
                // Synthetic cell VF sampled at tetra centers, not geometric intersection fractions.
                fractions.push_back(std::clamp(.5 - distance / step, 0., 1.));
            }
            for (const auto method : {TransferMethod::CellCenterSampling, TransferMethod::NodalInterpolation, TransferMethod::VolumeWeighted}) {
            // for (const auto method : {TransferMethod::CellCenterSampling, TransferMethod::NodalInterpolation}) {
                const std::string name = method == TransferMethod::CellCenterSampling ? "center" :
                    method == TransferMethod::NodalInterpolation ? "nodes" : "volumes";
                ExtractionOptions options;
                options.transferMethod = method;
                double extractionSeconds = 0;
                options.extractionSeconds = &extractionSeconds;
                const auto start = std::chrono::steady_clock::now();
                const Surface surface = extractVolumeSurface(points, offsets, nodes, fractions, {100,100,100}, options);
                const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                if (surface.triangles.empty()) throw std::runtime_error("Empty surface: " + shape + " / " + name);
                std::ofstream obj(output / (shape + "_" + name + ".obj"));
                obj.precision(17);
                for (const auto& vertex : surface.vertices) {
                    for (double coordinate : vertex)
                        if (!std::isfinite(coordinate)) throw std::runtime_error("Nonfinite vertex");
                    obj << "v " << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << '\n';
                }
                for (const auto& face : surface.triangles) {
                    for (auto index : face)
                        if (index >= surface.vertices.size()) throw std::runtime_error("Invalid surface index");
                    obj << "f " << face[0]+1 << ' ' << face[1]+1 << ' ' << face[2]+1 << '\n';
                }
                obj.flush();
                if (!obj) throw std::runtime_error("Cannot write OBJ");
                report << shape << ',' << name << ',' << surface.vertices.size() << ',' << surface.triangles.size()
                       << ',' << extractionSeconds << ',' << seconds << '\n';
                report.flush();
                if (!report) throw std::runtime_error("Cannot write report");
                std::cout << shape << " / " << name << ": " << surface.triangles.size() << " triangles, "
                          << seconds << " s" << std::endl;
            }
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
